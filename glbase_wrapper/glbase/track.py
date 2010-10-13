"""
track, part of glbase

TODO:
-----
. store the name (and other attribs?) in the db

"""

from __future__ import division

import cPickle, sys, os, struct, math, sqlite3, zlib

from location import location
from data import positive_strand_labels, negative_strand_labels

from numpy import array, zeros # Carray

TRACK_CACHE_SIZE = 10 # number of track segments to cache.

class track:
    """
    track definition, used for things like sequence reads across the genome
    """
    def __init__(self, name="None", new=False, filename=None, **kargs):
        """
        **Arguments**
            name (string)
                name for the track (defaults to filename)

            filename (string)
                directory location of the track file.
                only respected if dir_name is set.

        """
        if name:
            self.name = name
        else:
            self.name = filename

        # set-up the tables
        if new:
            self.__setup_tables(filename) # return an empty track
        else:
            self.__load_tables(filename)
        self.__c = None

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
        if path and not os.path.exists(path):
            os.makedirs(path)

        if os.path.exists(filename): # overwrite old file.
            os.remove(filename)
            # This could potentially fail - I should report and fail
            # nicely... At the moment it just throws an exception.

        self.__connection = sqlite3.connect(filename)
        self.__connection.text_factory = sqlite3.OptimizedUnicode

        c = self.__connection.cursor()

        c.execute("CREATE TABLE main (chromosome TEXT PRIMARY KEY, seq_reads INT)")

        self.__connection.commit()
        c.close()

    def __add_chromosome(self, chromosome):
        """
        add a chromosome to the main table.
        add a chromosome table.

        returns True if succesfully created/already present.
        """
        c = self.__connection.cursor()
        # check chromosome is not already present.
        if self.__has_chromosome(chromosome):
            return(True)

        c.execute("INSERT INTO main VALUES (?, ?)", (chromosome, 0)) # add chr to master table.

        # make the new chromsome table:
        table_name = "chr_%s" % str(chromosome)
        c.execute("CREATE TABLE %s (left INT, right INT, strand TEXT)" % (table_name, ))

        # link new table to old table.
        # how do I do that?!?!?

        c.close()
        return(True)

    def __has_chromosome(self, chromosome):
        """
        do we have that chromosome?
        """
        c = self.__connection.cursor()

        c.execute("SELECT chromosome FROM main WHERE chromosome=?", (chromosome, ))

        result = c.fetchone()

        c.close()

        if result:
            return(True)
        return(False)

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
        if not self.__has_chromosome(loc["chr"]):
            self.__add_chromosome(loc["chr"])

        c = self.__connection.cursor()

        # insert location into new array.
        table_name = "chr_%s" % str(loc["chr"])

        # get the old number of seq_reads
        c.execute("SELECT seq_reads FROM main WHERE chromosome=?", (loc["chr"], ))
        current_seq_reads = c.fetchone()[0] # always returns a tuple

        c.execute("UPDATE main SET seq_reads=? WHERE chromosome=?", (current_seq_reads+1, loc["chr"]))

        # add the location to the seq table:
        c.execute("INSERT INTO %s VALUES (?, ?, ?)" % table_name, (loc["left"], loc["right"], strand))

        c.close()

    def get_array(self, loc, strand=None, resolution=1, read_extend=0, **kargs):
        """
        **Purpose**
            get the data between location 'loc' and return it formatted as
            a nbp resolution array

        **Arguments**
            loc (Required)
                a valid location or string location.

            strand (Optional, default = None, ie. collect and merge both strands)
                strand, but only valid for stranded tracks
                if "+" return only that strand, if '-' return only the negative
                strand (will recognise several forms of strand, e.g. F/R, +/-

            resolution (Optional, default = 1bp)
                nbp resolution required (you should probably send a float for accurate rendering)

            read_extend (Optional, default = 0)
                extend the read length to 'fill in the peak'
                if the original reads are 36bp, then add ~70bp to give an
                estimated size of the peak.
                If the reads are end-based, then set this to the estimated
                size of the DNA shear.

        **Returns**
            an 'numpy.array([0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
        # check if the array is already on the cache.
        extended_loc = location(loc=str(loc)) # this should be fixed in a later version of glbase...
        extended_loc.expand(read_extend) # this behaviour may not work in a future version?
        locs_required = [extended_loc] # make sure to include the read_extend, but don't modify the location to get.

        # get the reads only in the ranges required.
        # at the moment this is not implemented.
        reads = []
        for l in locs_required:
            reads += self.get_reads(l, strand)

        # make a single array
        a = zeros(int( (loc["right"]-loc["left"]+resolution)/resolution ))

        # work out and standardise strand:
        if strand:
            if strand in positive_strand_labels:
                strand = positive_strand_labels
            elif strand in negative_strand_labels:
                strand = negative_strand_labels

        for r in reads:
            if not strand or strand and r[2] in strand:
                if r[2] in positive_strand_labels:
                    read_left = r[0]
                    read_right = r[1] + read_extend
                elif r[2] in negative_strand_labels:
                    read_left = r[0] - read_extend
                    read_right = r[1]

                for rloc in xrange(read_left, read_right, int(resolution)):
                    array_relative_location = int((rloc - read_extend - loc["left"]) / resolution) # convert relative to the array

                    if array_relative_location >= 0 and array_relative_location < len(a): # within array
                        a[array_relative_location] += 1

        return(a)

    def get_array_chromosome(self, chromosome, strand=None, resolution=1, read_extend=0, **kargs):
        """
        **Purpose**
            get the entire array data for the chromosome

        **Arguments**
            chromosome (Required)
                a number '1', string "1", or X, Y

            strand (Optional, default = None, ie. collect and merge both strands)
                strand, but only valid for stranded tracks
                if "+" return only that strand, if '-' return only the negative
                strand (will recognise several forms of strand, e.g. F/R, +/-

            resolution (default = 1bp)
                nbp resolution required (you should probably send a float for accurate rendering)

            read_extend (Optional, default = 0)
                extend the read length to 'fill in the peak'
                if the original reads are 36bp, then add ~70bp to give an
                estimated size of the peak.
                If the reads are end-based, then set this to the estimated
                size of the DNA shear.

        **Returns**
            an 'array('i', [0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
        if strand: raise NotImplementedError, "Eh... strand not supported yet..."

        if not self.__c:
            self.__c = self.__connection.cursor()
        table_name = "chr_%s" % chromosome
        self.__c.execute("SELECT * FROM %s" % table_name)
        reads = self.__c.fetchall()

        # I need to find the right most read to estimate the size of the track array.
        right_most = 0
        for i in reads:
            if right_most < i[1]+read_extend:
                right_most = i[1]+read_extend

        # make an array.
        a = zeros(int(right_most+int(resolution), int(resolution)))

        for r in reads:
            # read_extend
            if r[2] in positive_strand_labels:
                left = r[0]
                right = r[1] + read_extend
            elif r[2] in negative_strand_labels:
                left = r[0] - read_extend
                right = r[1]

            # check for array wrap arounds:
            if left < 0: left = 0

            for loc in xrange(left, right, 1):
                if resolution <= 1: # force 1bp read # this may be incorrect as it may add the read several times until it increments?
                    # this is the source of the artivacts when resolution < 1.0?
                    a[loc] += 1
                else:
                    a[int(loc/resolution)] += 1

        #print "array_len", len(a)

        return(a)

    def get_reads(self, loc, strand=None):
        """
        **Purpose**
            get all of the sequence reads between location 'loc' and return
            it formatted as a list of tuples: (left, right, strand), seq reads.

        **Arguments**
            loc (Required)
                a valid location or string location.

        **Returns**
            a list containing all of the reads between loc.
        """
        if not self.__c:
            self.__c = self.__connection.cursor()

        table_name = "chr_%s" % loc["chr"]

        # ~4.11 s 946 reads
        self.__c.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
            (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))
        #self.__c.execute("SELECT left, right, strand FROM %s WHERE (left<=? AND right>=?) OR (?<=left AND ?>=right) OR (?>=left AND ?<=right) OR (?>=left AND ?<=right)" % table_name,
        #    (loc["left"], loc["right"], loc["left"], loc["right"], loc["left"], loc["left"],loc["right"], loc["right"]))
        # pseudo code:
        # if left - dbright < 0 and :
        #span = loc["right"] - loc["left"]
        #self.__c.execute("SELECT left, right, strand FROM %s WHERE (?-left > -1 AND ?-left<?) OR (?-right>-1 AND ?-right<?)" % table_name,
        #    (loc["right"], loc["right"], span, loc["left"], loc["left"], span))

        result = self.__c.fetchall()

        return(result)

    def get_read_count(self, loc):
        """
        **Purpose**
            get the number of reads within the location specified

        **Arguments**
            loc (Required)
                a valid location or string location.

        **Returns**
            an integer (or 0) containing the number of reads falling within
            the location string.
        """
        if not self.__c:
            self.__c = self.__connection.cursor()

        table_name = "chr_%s" % loc["chr"]

        self.__c.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
            (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))

        # this is here in the hope that in the future cursor.rowcount
        # will be correctly supported...
        # at the moment probably provides little or no benefit.

        return(len(self.__c.fetchall()))

    def finalise(self):
        """
        finalise the database (shrink unused edit space)
        dump useless bits etc.
        You must call this! to finalise the db.
        Or get() will not work!
        This copies the cache onto disk and closes the db.
        """
        # do a commit
        self.__connection.commit()
        self.__connection.execute("VACUUM")
        self.__connection.commit()

    def _debug__print_all_tables(self):
        c = self.__connection.cursor()
        c.execute("SELECT * FROM main")
        result = c.fetchall()

        print "Main:"
        for item in result:
            print item

        print "Chr_Tables:"
        for item in result:
            table_name = "chr_%s" % str(item[0])[0] # stop injections.
            print " Table", table_name
            c.execute("SELECT * FROM %s" % table_name)
            chr_table_res = c.fetchall()
            for i in chr_table_res:
                print " ", i
        c.close()

if __name__ == "__main__":
    # a few testers to see how it works.

    """
    t = track(filename="data/SpMash1_new.trk")
    # speed tester

    from random import randint
    import time
    from numpy import mean, std

    chroms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"X", "Y"]
    results = []
    lens = [] # sanity checking.

    for s in xrange(10): # number of samples
        s = time.time()
        for n in xrange(10):
            r = t.get_reads(location(loc="chr5:114423935-114444335"))# % (chroms[randint(0, len(chroms)-1)])))
        e = time.time()
        results.append(e - s)
    print "Took: mean %s +- %s" % (mean(results), std(results))
    print "number of reads:", len(r)

    """
    t = track(filename="data/test_new.trk", new=True)
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

    print "Finalise:"
    t.finalise() # must call this
    t._debug__print_all_tables()

    # see the test_case for a formal test.
    print t.get_reads(location(loc="chr1:1-26"))
    print t.get_array(location(loc="chr1:1-26"))
    print t.get_array(location(loc="chr1:10-20"))
    print t.get_array(location(loc="chr1:20-25"))

    print "\nReload"
    t = track(filename="data/test_new.trk")
    print t.get_reads(location(loc="chr1:1-26"))
    print t.get_array(location(loc="chr1:1-26"))
    print t.get_array(location(loc="chr1:10-20"))
    print "res_tests"
    print t.get_array(location(loc="chr1:20-30"))
    print t.get_array(location(loc="chr1:20-30"), resolution=2)
    print t.get_array(location(loc="chr1:20-30"), resolution=1.5)
    #print t.get_array(location(loc="chr1:20-30"), resolution=0.5) # this causes an error
    print t.get_array(location(loc="chr1:20-30"), resolution=5)

    #print t.get_array_chromosome("1")
