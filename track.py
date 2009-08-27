"""
track, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

"""

from __future__ import division

import cPickle, sys, os, struct, ConfigParser, math, sqlite3, zlib # renamed to configparser in >2.6

from glbase_wrapper import location as location
from glbase_wrapper import positive_strand_labels, negative_strand_labels

from array import array

TRACK_BLOCK_SIZE = 1000 # should go in opt, later

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

        self.__connection = sqlite3.connect(filename)
        c = self.__connection.cursor()

        c.execute("CREATE TABLE data (blockID TEXT PRIMARY KEY, chrom INTEGER, left INTEGER, right INTEGER, array BLOB)")

        self.__connection.commit()
        # set up the block table describing the number of blocks present.
        c.execute("CREATE TABLE info (chrom INTEGER PRIMARY KEY, blockSize INTEGER)")

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
            True
        """

        left_most_block = int(abs(math.floor(loc["left"] / self.block_size)))
        right_most_block = int(abs(math.ceil((loc["right"]+1) / self.block_size)))

        blocks_required = ["%s:%s" % (loc["chr"], b) for b in xrange(left_most_block * self.block_size, right_most_block * self.block_size, self.block_size)]

        for blockID in blocks_required:
            # this is the span location of the block
            block_loc = location(chr=blockID.split(":")[0], left=blockID.split(":")[1], right=int(blockID.split(":")[1])+self.block_size-1)

            # check if the block already exists
            if not self.__has_block(blockID): # not present, make a new one.
                self.__new_block(block_loc["chr"], block_loc["left"])
            block = self.__get_block(blockID)

            array_data = self.__unformat_data(block[1]) # get back the usable data

            # modify the data
            for pos in xrange(self.block_size): # iterate through the array.
                local_pos = block_loc["left"] + pos
                # some funny stuff here:
                # right is inc'd by 1
                # as that is what is usually expected from 10-20.
                # (i.e. coords are NOT 0-based and are closed).
                if local_pos < (loc["right"]+1): # still within block
                    if local_pos >= loc["left"] and local_pos <= (loc["right"]+1): # within the span to increment.
                        array_data[pos] += increment

            self.__modify_block(block, array_data) # load it back in.
        #print "\n\nDebug:"
        #self.__see_block_counts()
        #self.__see_data()
        return(True)

    def __has_block(self, blockID):
        """
        checks if the data base has that block already.
        returns only True or False,
        does not return the block, you must use __get_block()
        to get the actual block.
        """
        has = False
        c = self.__connection.cursor()
        c.execute("SELECT blockID FROM data WHERE blockID=?", (blockID, ))
        result = c.fetchone()
        if result:
            has = True
        c.close()
        return(has)

    def __modify_block(self, old_block, new_data):
        """
        update the block with new data.

        block = [blockID, chr, left, right, data]
        """
        # get the block data.
        c = self.__connection.cursor()
        c.execute("UPDATE data SET array=? WHERE blockID=?", (self.__format_data(new_data), old_block[0]))
        self.__connection.commit()
        c.close()

    def __new_block(self, chr, left, data=None):
        """
        add a data block to the db in data table.
        returns only True, does not return the block.
        use self.__get_block() to get a block.
        """
        if not data: # fill a blank entry
            data = array('i', [])
            for n in xrange(self.block_size):
                data.append(0)

        c = self.__connection.cursor()

        blockID = "%s:%s" % (chr, left)

        c.execute("INSERT INTO data VALUES (?,?,?,?,?)", (blockID, chr, left, left+self.block_size-1, self.__format_data(data)))

        # get old value:
        c.execute("SELECT blockSize FROM info WHERE chrom=?", (chr, ))
        result = c.fetchone()

        if result: # already present, modify:
            c.execute("UPDATE info SET blockSize=? WHERE chrom=?", (result[0]+1, chr)) # result[0] = old block size
        else:
            c.execute("INSERT INTO info VALUES (?,?)", (chr, 1)) # make a new block

        self.__connection.commit()
        c.close()

        return(True)

    def __get_block(self, blockID):
        """
        get the block identified by chr and left coordinates and return a Python Object.
        """
        c = self.__connection.cursor()
        c.execute("SELECT blockID, array FROM data WHERE blockID=?", (blockID, ))
        result = c.fetchone()
        c.close()
        if result:
            return(result)
        else:
            raise Exception, "No Block!"

    def __format_data(self, data):
        """
        array('i', []) --> whatever it's stored as in db
        """
        return(str(data))

    def __unformat_data(self, data):
        """
        whatever stored as in db --> array('i', [])
        """
        return(eval(data))

    def __see_block_counts(self):
        c = self.__connection.cursor()
        c.execute("SELECT * FROM info")
        print "Info contents:-----"
        for line in c:
            print line
        print "End----------------"
        c.close()

    def __see_data(self):
        c = self.__connection.cursor()
        c.execute("SELECT * FROM data")
        print "Data contents:-----"
        for line in c:
            print line
        print "End----------------"
        c.close()

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

        blocks_required = ["%s:%s" % (loc["chr"], b) for b in xrange(left_most_block * self.block_size, right_most_block * self.block_size, self.block_size)]

        ret_array = array('i', [])

        for blockID in blocks_required:
            # this is the span location of the block
            block_loc = location(chr=blockID.split(":")[0], left=blockID.split(":")[1], right=int(blockID.split(":")[1])+self.block_size-1)

            # check if the block already exists
            if self.__has_block(blockID): # it does, get it.
                block = self.__get_block(blockID)
                this_block_array_data = self.__unformat_data(block[1]) # get back the usable data
            else: # not present, fake a block.
                this_block_array_data = array('i', [0 for x in xrange(self.block_size)])

            # modify the data
            for pos in xrange(self.block_size): # iterate through the array.
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
        sort the list by Chrom > left.
        """
        c = self.__connection.cursor()
        c.execute("VACUUM")
        self.__connection.commit()

    def show_debug_info(self):
        self.__see_block_counts()

if __name__ == "__main__":
    # a few testers to see how it works.

    t = track(filename="data/test.trk", new=True)
    #t.add_location(location(loc="chr1:10-20"))
    t.add_location(location(loc="chr1:10-22"))
    t.add_location(location(loc="chr1:10-21"))
    t.add_location(location(loc="chr1:10-20"))
    t.add_location(location(loc="chr1:10-19"))
    #t.add_location(location(loc="chr1:1-31"))

    #print t.get("chr1:1-35")

    t.add_location(location(loc="chr1:11-21"))
    t.add_location(location(loc="chr1:12-22"))
    t.add_location(location(loc="chr1:25-26"))
    t.add_location(location(loc="chr1:25-26"))
    t.add_location(location(loc="chr1:25-26"))
    # really boost one single location:

    for n in xrange(50):
        t.add_location(location(loc="chr1:11-12"))
        t.add_location(location(loc="chr1:12-13"))

    t.add_location(location(loc="chr1:10-20"), strand="-")
    t.add_location(location(loc="chr1:11-21"), strand="-")
    t.add_location(location(loc="chr1:12-22"), strand="-")
    t.add_location(location(loc="chr1:25-26"), strand="-")

    t.add_location(location(loc="chr2:10-20"), strand="-")
    t.add_location(location(loc="chr2:11-21"), strand="-")
    t.add_location(location(loc="chr2:12-22"), strand="-")
    t.add_location(location(loc="chr2:25-26"), strand="-")

    print t.get(location(loc="chr1:1-26"))
    print t.get(location(loc="chr1:10-20"))
    print t.get(location(loc="chr1:20-25"))
    print t.get("chr1:10-20")
    t.finalise()

    print "\nReload"
    t = track(filename="data/test.trk")
    print t.get(location(loc="chr1:1-26"))
    print t.get(location(loc="chr1:10-20"))
    print t.get(location(loc="chr1:20-25"))
    print t.get("chr1:10-20")
