"""
track, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

"""



import pickle, sys, os, struct, configparser, math, sqlite3, zlib # renamed to configparser in >2.6

from glbase_wrapper import location as location
from glbase_wrapper import positive_strand_labels, negative_strand_labels

from array import array

TRACK_BLOCK_SIZE = 2000 # should go in opt, later
CACHE_SIZE = 100000 # maximum number of blocks to keep in memory.
# these are equivalent to about 800 Mb's of memory on a Fedora box
# On a 2Gb machine it is still usable, less may be a problem.
# needs to be tested on a Windows and MACOSX machine.

class track:
    """
    track definition, used for things like sequence reads across the genome
    """
    def __init__(self, name="None", new=False, stranded=True, filename=None, bin_format="i"):
        """
        **Arguments**
            name (string)
                name for the track (defaults to "None")

            stranded (True|False)
                store one or both strands?

            filename (string)
                directory location of the track file.
                only respected if dir_name is set.

        """
        self.name = name
        self.bin_format = bin_format

        bin_len = {"i": struct.calcsize("i"), "l": struct.calcsize("l")}
        self.bin_len = bin_len[self.bin_format]

        self.block_size = TRACK_BLOCK_SIZE # size of the blocks

        self.strands = ["+"]
        if stranded: self.strands += ["-"]

        self.cache = {}
        self.cacheQ = []

        # set-up the tables
        if new:
            self.__setup_tables(filename) # return an empty track
        else:
            self.__load_tables(filename)

    def __load_tables(self, filename):
        """
        just load in the tables.
        (basically, fill __connection)
        """
        self.__connection = sqlite3.connect(filename)

    def __setup_tables(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # kill any previously exisiting file (Use with care!)
        # make sure the directory is available:
        path = os.path.split(filename)[0]
        if not os.path.exists(path):
            os.makedirs(path)

        if os.path.exists(filename): # overwrite old file.
            os.remove(filename)
            # This could potentially fail - I should report and fail
            # nicely... At the moment it just throws an exception.

        self.__connection = sqlite3.connect(filename)
        self.__connection.text_factory = sqlite3.OptimizedUnicode

        c = self.__connection.cursor()

        c.execute("CREATE TABLE data (blockID TEXT PRIMARY KEY, array TEXT)")

        self.__connection.commit()
        c.close()

    def add_location(self, loc, strand="+", increment=1):
        """
        **Purpose**
            Add a location to the track.
            Increments the score by 'increment' from loc["left"] to
            loc["right"]

        **Arguments**
            loc

            strand

            increment

        **Returns**
            True, if completes succesfully, or exception.
        """
        left_most_block = int(abs(math.floor(loc["left"] / self.block_size)))
        right_most_block = int(abs(math.ceil((loc["right"]+1) / self.block_size)))

        blocks_required = ["%s:%s" % (loc["chr"], b) for b in range(left_most_block * self.block_size, right_most_block * self.block_size, self.block_size)]

        for blockID in blocks_required:
            # this is the span location of the block
            #block_loc = location(chr=blockID.split(":")[0], left=blockID.split(":")[1], right=int(blockID.split(":")[1])+self.block_size-1)

            # check if the block already exists
            if not self.__has_block(blockID): # not present, make a new one.
                self.__new_block(blockID)
            else:
                if blockID not in self.cacheQ: # not on cache, add it;
                    self.cacheQ.insert(0, blockID) # put the ID at the front.
                    self.cache[blockID] = self.__get_block(blockID)
            # block should now be on cache and accesable.

            bleft = int(blockID.split(":")[1])
            lleft = int(loc["left"])
            lright = int(loc["right"])
            # modify the data
            for pos in range(self.block_size): # iterate through the array.
                local_pos = bleft + pos # only require "left"
                # some funny stuff here:
                # right is inc'd by 1
                # as that is what is usually expected from 10-20.
                # (i.e. coords are NOT 0-based and are closed).
                if local_pos >= lleft and local_pos <= lright: # within the span to increment.
                    self.cache[blockID][pos] += increment

            self.__flush_cache()
        return(True)

    def __has_block(self, blockID):
        """
        checks if the data base has that block already.
        returns only True or False,
        does not return the block, you must use __get_block()
        to get the actual block.
        """
        if blockID in self.cache:
            return(True) # on cache, so must exist.

        has = False
        c = self.__connection.cursor()
        c.execute("SELECT blockID FROM data WHERE blockID=?", (blockID, ))
        result = c.fetchone()
        if result:
            has = True
        c.close()
        return(has)

    def __commit_block(self, blockID, data):
        """
        update the block with new data.
        commit block to db.
        """
        # see if tyhe block is in the database:
        c = self.__connection.cursor()
        c.execute("SELECT blockID FROM data WHERE blockID=?", (blockID, ))
        result = c.fetchone()

        d = self.__format_data(data)

        if result: # has a block already, modify it.
            # update the block data.

            c.execute("UPDATE data SET array=? WHERE blockID=?", (self.__format_data(data), blockID))
        else:
            c.execute("INSERT INTO data VALUES (?, ?)", (blockID, d))
        c.close()

    def __new_block(self, blockID, data=None):
        """
        add a data block to the db in data table.
        returns only True, does not return the block.
        use self.__get_block() to get a block.

        new_block DOES NOT WRITE INTO THE DB!
        You need to flush the cache for that to happen
        """
        if not data: # fill a blank entry
            data = array('i', [0 for x in range(self.block_size)])
            #data.extend([1 for x in xrange(self.block_size)])
            #for n in xrange(self.block_size):
            #    data.append(0)

        if blockID not in self.cacheQ: # not on cache
            self.cacheQ.insert(0, blockID) # put the ID at the front.
        self.cache[blockID] = data

        return(False)

    def __flush_cache(self):
        """
        check the cache is not over the size limit. If it is, take the last
        n>CACHE_SIZE entries commit to the db.

        """
        while len(self.cacheQ) > CACHE_SIZE:
            blockID = self.cacheQ.pop()
            self.__commit_block(blockID, self.cache[blockID])
            del self.cache[blockID]

        return(True)

    def __get_block(self, blockID):
        """
        get the block identified by chr and left coordinates and return a Python Object.
        """
        if blockID in self.cache:
            return(self.cache[blockID])

        # not on the cache. get the block and put it on the cache.
        c = self.__connection.cursor()
        c.execute("SELECT array FROM data WHERE blockID=?", (blockID, ))
        result = c.fetchone()
        c.close()

        if result:
            return(self.__unformat_data(result[0]))
        else:
            raise Exception("No Block!")

    def __format_data(self, data):
        """
        array('i', []) --> whatever it's stored as in db
        """
        return(sqlite3.Binary(zlib.compress(data.tostring())))

    def __unformat_data(self, data):
        """
        whatever stored as in db --> array('i', [])
        """
        #print "ret:",[d for d in data], ":"
        a = array('i')
        a.fromstring(zlib.decompress(data))
        return(a)

    def get(self, loc, strand="+"):
        """
        **Purpose**
            get the data between location 'loc'

        **Arguments**
            loc (Required)
                a valid location or string location.

            strand (Optional, default = "+")
                strand, but only valid for stranded tracks

        **Returns**
            an 'array('i', [0, 1, 2 ... n])' contiginous array
        """
        try:
            if loc["chr"]: pass
        except TypeError: # probably a location string. try to cooerce
            loc = location(loc=loc)

        left_most_block = int(abs(math.floor(loc["left"] / self.block_size)))
        right_most_block = int(abs(math.ceil((loc["right"]+1) / self.block_size)))

        blocks_required = ["%s:%s" % (loc["chr"], b) for b in range(left_most_block * self.block_size, right_most_block * self.block_size, self.block_size)]

        ret_array = array('i', [])

        for blockID in blocks_required:
            # this is the span location of the block
            block_loc = location(chr=blockID.split(":")[0], left=blockID.split(":")[1], right=int(blockID.split(":")[1])+self.block_size-1)

            # check if the block already exists
            if self.__has_block(blockID): # it does, get it.
                block = self.__get_block(blockID)
                this_block_array_data = block # get back the usable data
            else: # block not in db, fake a block instead.
                this_block_array_data = array('i', [0 for x in range(self.block_size)])

            #print "b", block
            #print self.__get_block(blockID)

            # modify the data
            for pos in range(self.block_size): # iterate through the array.
                local_pos = block_loc["left"] + pos
                # see add_location for details
                if local_pos < (loc["right"]+1): # still within block
                    if local_pos >= loc["left"] and local_pos <= (loc["right"]+1): # within the span to increment.
                        ret_array.append(this_block_array_data[pos])
        return(ret_array)

    def finalise(self):
        """
        finalise the database (shrink unused edit space)
        dump useless bits etc.
        You must call this! to finalise the db.
        Or get() will not work!
        This copies the cache onto disk and closes the db.
        """
        for blockID in self.cache:
            self.__commit_block(blockID, self.cache[blockID])
        self.cache = {}
        self.cacheQ = []
        self.__connection.commit() # commit all the __commit_block()

        c = self.__connection.cursor()
        c.execute("VACUUM")
        self.__connection.commit()

if __name__ == "__main__":
    # a few testers to see how it works.

    t = track(filename="data/test.trk", new=True)
    t.add_location(location(loc="chr1:1-30"))
    t.add_location(location(loc="chr1:1-30"))
    t.add_location(location(loc="chr1:10-22"))
    t.add_location(location(loc="chr1:10-21"))
    t.add_location(location(loc="chr1:10-20"))
    t.add_location(location(loc="chr1:10-19"))
    t.add_location(location(loc="chr1:11-21"))
    t.add_location(location(loc="chr1:12-22"))
    t.add_location(location(loc="chr1:24-25"))
    t.add_location(location(loc="chr1:25-26"))
    t.add_location(location(loc="chr1:26-27"))
    t.add_location(location(loc="chr1:26-26"))
    t.add_location(location(loc="chr1:27-27"))

    t.add_location(location(loc="chr1:1-100"), strand="+")
    t.add_location(location(loc="chr1:26-100"), strand="+")

    t.add_location(location(loc="chr2:100-200"), strand="-")
    t.add_location(location(loc="chr2:110-210"), strand="-")
    t.add_location(location(loc="chr2:120-220"), strand="-")
    t.add_location(location(loc="chr2:250-260"), strand="-")

    print("Finalise:")
    t.finalise() # must call this

    print(t.get(location(loc="chr1:1-26")))
    print(t.get(location(loc="chr1:10-20")))
    print(t.get(location(loc="chr1:20-25")))

    print("\nReload")
    t = track(filename="data/test.trk")
    print(t.get(location(loc="chr1:1-26")))
    print(t.get(location(loc="chr1:10-20")))
    print(t.get(location(loc="chr1:20-25")))
