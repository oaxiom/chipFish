"""
track, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

TODO:
-----
. store the name (and other attribs?) in the db

"""

from __future__ import division

import cPickle, sys, os, struct, math, sqlite3, zlib

from glbase_wrapper import location as location
from glbase_wrapper import positive_strand_labels, negative_strand_labels

from array import array

class track:
    """
    track definition, used for things like sequence reads across the genome
    """
    def __init__(self, name="None", new=False, stranded=True, filename=None, array_format="i"):
        """
        **Arguments**
            name (string)
                name for the track (defaults to filename)

            stranded (True|False)
                store one or both strands?

            filename (string)
                directory location of the track file.
                only respected if dir_name is set.

        """
        if name:
            self.name = name
        else:
            self.name = filename

        self.array_format = array_format

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

            resolution (Not Implemented)
                nbp resolution required (you should probably send a float)

            read_extend (Optional, default = 0, Not Implemented)

        **Returns**
            an 'array('i', [0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
        # check if the array is already on the cache.
        locs_required = [loc]

        # get the reads only in the ranges required.
        reads = []
        for l in locs_required:
            reads += self.get_reads(l, strand)

        #print resolution
        #print len(xrange(0, loc["right"]-loc["left"]+int(resolution), resolution))

        # make an array
        a = array(self.array_format, [0 for x in xrange(0, loc["right"]-loc["left"]+int(resolution), int(resolution))])

        for loc in locs_required:
            for real_loc in xrange(loc["left"], int(loc["right"]+resolution), int(resolution)):
                for r in reads:
                    if r[2] in positive_strand_labels:
                        left = r[0]
                        right = r[1] + read_extend
                    elif r[2] in negative_strand_labels:
                        left = r[0] - read_extend
                        right = r[1]

                    if real_loc >= left and real_loc <= right:
                        a[int((real_loc-loc["left"])/resolution)] += 1

        #print "array_len", len(a)

        return(a) # py3k problem here...

    def get_reads(self, loc, strand=None):
        """
        **Purpose**
            get all of the sequence reads between location 'loc' and return
            it formatted as a list of tuples: (left, right, strand), seq reads.

        **Arguments**
            loc (Required)
                a valid location or string location.

            strand (Optional, default = 'None', ie. collect both strands)
                collect only one strands data. specify + or - strand and
                return only seq reads on that strand.

        **Returns**
            a list containing all of the reads between loc.
        """
        c = self.__connection.cursor()

        table_name = "chr_%s" % loc["chr"]

        if not strand:
            c.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
                (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))
        else:
            pass

        result = c.fetchall()

        c.close()
        return(result)

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
    print t.get_array(location(loc="chr1:20-25"))
    print t.get_array(location(loc="chr1:20-25"), resolution=2)
    print t.get_array(location(loc="chr1:20-25"), resolution=1.5)

    print t.get_array(location(loc="chr1:20-25"), resolution=1.5)
