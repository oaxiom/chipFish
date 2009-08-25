"""
track, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

"""

from __future__ import division

import cPickle, sys, os, struct, ConfigParser, math, sqlite3 # renamed to configparser in >2.6

from glbase_wrapper import location as location
from glbase_wrapper import positive_strand_labels, negative_strand_labels

from array import array

TRACK_BLOCK_SIZE = 10 # should go in opt, later

class track:
    """
    track definition, used for things like sequence reads across the genome
    """
    def __init__(self, name="None", new=True, stranded=True, filename=None, bin_format="i"):
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
        self.track = {}
        self.bin_format = bin_format

        bin_len = {"i": struct.calcsize("i"), "l": struct.calcsize("l")}
        self.bin_len = bin_len[self.bin_format]

        self.block_size = TRACK_BLOCK_SIZE # size of the blocks

        self.strands = ["+"]
        if stranded: self.strands += ["-"]

        # set-up the tables
        if new:
            self._setup_tables(filename) # return an empty track
        else:
            self._load_tables(filename)

    def _setup_tables(self, filename):
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

        array_data = array('i', []) # build a block
        for item in xrange(self.block_size):
            array_data.append(0)
        c.execute("INSERT INTO data VALUES (?,?,?,?,?)", ["1:11", 1, 11, 20, str(array_data)])

        c.execute("INSERT INTO data VALUES (?,?,?,?,?)", ["1:21", 1, 21, 30, ""])
        c.execute("INSERT INTO data VALUES (?,?,?,?,?)", ["1:31", 1, 31, 40, str(array_data)])

        self.__connection.commit()

        c.execute("SELECT * FROM data")

        for item in c:
            print ":", item

        print "\n>block1"
        c.execute("SELECT * FROM data WHERE blockID='1:11'")

        for item in c:
            print item

        print "\n>block2"
        c.execute("SELECT * FROM data WHERE blockID='1:21'")

        for item in c:
            print item

        c.close()

        # set up the block table describing the number of blocks present.
        c.execute("CREATE TABLE blocksize (chrom INTEGER PRIMARY KEY, blockSize INTEGER)")


    def add_location(self, loc, strand="+"):
        c = self.__connection.cursor()

        left_most_block = abs(math.floor(loc["left"] / self.block_size)
        right_most_block= abs(math.ceil(loc["right"] / self.block_size)
        
        blocks_required = ["%s:%s" (loc["chr"], b) for b in xrange(left_most_block, right_most_block)]

        array_data = array('i', []) # build a block
        for item in xrange(self.block_size):
            array_data.append(0)
        c.execute("INSERT INTO data VALUES (?,?,?,?,?)", ["%s:%s" % (loc["chr"], loc["left"]), 1, 11, 20, str(array_data)])

        self.__connection.commit()

    def __getitem__(self, value):
        pass

if __name__ == "__main__":
    # a few testers to see how it works.

    t = track(filename="data/test/out.trk", new=True)
    t.add_location(location(loc="chr1:10-20"))
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

    print t.track
    print t[location(loc="chr1:1-26")]
    print t[location(loc="chr1:10-20")]
    print t[location(loc="chr1:20-25")]
    print t["chr1:10-20"]

    #print t[1] # error

    t = track(delayed=True, dir_name="data/test/")
    print t[location(loc="chr1:1-26")]
    print t[location(loc="chr1:10-20")]
    print t[location(loc="chr1:20-25")]
    print t["chr1:10-20"]
