"""

An SQL replacement for genome

This loses a lot of the flexibility of 

"""

from __future__ import division

import cPickle, sys, os, struct, math, sqlite3, zlib, time, csv

from operator import itemgetter

from progress import progressbar
from errors import AssertionError
from location import location
import utils, config
from data import positive_strand_labels, negative_strand_labels
from draw import draw
import matplotlib.cm as cm
import matplotlib.pyplot as plot
import scipy.stats as stats
from scipy.stats import pearsonr
from base_track import base_track

import numpy
from numpy import array, zeros, set_printoptions, int32, append, linspace, argmax, amax, delete

TRACK_CACHE_SIZE = 10 # number of track segments to cache.

class genome_sql(base_track):
    """
    track definition, used for things like sequence reads across the genome
    
    **Arguments**
        name (string)
            name for the track (defaults to filename)

        filename (string)
            directory location of the track file.
            only respected if dir_name is set.
        
        new (Optional, default=False)
            Use seqToTrk() in preference of this. But if you know what you are 
            doing then this will generate a new (empty) db.

    """
    def __init__(self, name=None, new=False, filename=None, norm_factor=1.0, **kargs):
        base_track.__init__(self, name, new, filename, norm_factor)
        
        if new:
            self.__setup_tables(filename)

    def __repr__(self):
        return("glbase.genome_sql")

    def __setup_tables(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # If i make it here, then a lot of grunt work already done in base_track
        c = self._connection.cursor()

        c.execute("CREATE TABLE main (chromosome TEXT PRIMARY KEY, num_features INT)")

        self._connection.commit()
        c.close()

    def __add_chromosome(self, chromosome):
        """
        add a chromosome to the main table.
        add a chromosome table.

        returns True if succesfully created/already present.
        """
        c = self._connection.cursor()
        # check chromosome is not already present.
        if self.__has_chromosome(chromosome):
            return(True)

        c.execute("INSERT INTO main VALUES (?, ?)", (chromosome, 0)) # add chr to master table.

        # make the new chromsome table:
        table_name = "chr_%s" % str(chromosome)
        c.execute("CREATE TABLE %s (transcript_left INT, transcript_right INT, cds_left INT, cds_right INT, exonStarts TEXT, exonEnds TEXT, name TEXT, strand TEXT)" % (table_name, ))

        c.close()
        return(True)

    def __has_chromosome(self, chromosome):
        """
        do we have that chromosome?
        """
        c = self._connection.cursor()

        c.execute("SELECT chromosome FROM main WHERE chromosome=?", (chromosome, ))

        result = c.fetchone()

        c.close()

        if result:
            return(True)
        return(False)

    def add_feature(self, loc, cds_loc, exonCounts, exonStarts, exonEnds, name, strand="+", increment=1):
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

        # convert strand to +, -
        if strand in positive_strand_labels:
            strand = "+"
        elif strand in negative_strand_labels:
            strand = "-"

        c = self._connection.cursor()

        # insert location into new array.
        table_name = "chr_%s" % str(loc["chr"])

        # get the old number of seq_reads
        c.execute("SELECT seq_reads FROM main WHERE chromosome=?", (loc["chr"], ))
        current_seq_reads = c.fetchone()[0] # always returns a tuple

        c.execute("UPDATE main SET seq_reads=? WHERE chromosome=?", (current_seq_reads+1, loc["chr"]))

        # add the location to the seq table:
        insert_data = (loc['left'], loc['right'], cds_loc['left'], cds_loc['right'], str(exonStarts), str(exonEnds), str(name), str(strand))
        c.execute("INSERT INTO %s VALUES (?, ?, ?)" % table_name, insert_data)

        c.close()

    def get(self, loc, resolution=1, read_extend=0, kde_smooth=False, 
        view_wid=0, strand=False, **kargs):
        """
        **Purpose**
            get the data between location 'loc' and return it formatted as
            a nbp resolution array

        **Arguments**
            loc (Required)
                a valid location or string location.

            resolution (Optional, default = 1bp)
                nbp resolution required (you should probably send a float for accurate rendering)

            read_extend (Optional, default = 0)
                extend the read length to 'fill in the peak'
                if the original reads are 36bp, then add ~70bp to give an
                estimated size of the peak.
                If the reads are end-based, then set this to the estimated
                size of the DNA shear.
            
            kde_smooth (Experimental)
                perform kde smoothng on the data, using the integer specified as an option.
                In this case the read_extend acts as a tag shift instead of a read_extend
                Hence set that to half of the expected shear size.

            strand (Optional, default=False)
                collect only reads on the specified strand. (track will use read strand 
                information intelligently, if present).
    
        **Returns**
            an 'numpy.array([0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
        if not isinstance(loc, location):
            loc = location(loc=loc)
        extended_loc = loc.expand(read_extend)

        result = self.get_reads(extended_loc, strand=strand)
        
        if kde_smooth:
            return(self.__kde_smooth(loc, result, resolution, 0, view_wid, read_extend))

        loc_left = loc["left"]
        loc_right = loc["right"]

        # make a single array
        a = [0] * int( (loc_right-loc_left+resolution)/resolution ) # Fast list allocation
        # Python lists are much faster for this than numpy or array

        len_a = len(a)
        
        for r in result:
            read_left, read_right, strand = r
            if strand == "+":
                read_right += (read_extend + 1) # coords are open
            elif strand == "-" :
                read_left -= read_extend
                read_right += 1 # coords are open 
            
            rel_array_left = int((read_left - loc_left) // resolution)
            rel_array_right = int((read_right - loc_left) // resolution)        
            
            if rel_array_left <= 0:
                rel_array_left = 0
            if rel_array_right > len_a:
                rel_array_right = len_a
            
            for array_relative_location in xrange(rel_array_left, rel_array_right, 1):
                a[array_relative_location] += 1
            
            #a[rel_array_left:rel_array_right] += 1 # Why is this slower than the for loop? # could be done with num_expr?
            
            #[a[array_relative_location].__add__(1) for array_relative_location in xrange(rel_array_left, rel_array_right, 1)] # just returns the exact item, a is unaffected?
        return(numpy.array(a)*self.norm_factor)             

    def getFeatures(self, loc=None, **kargs):
        """
        **Purpose**
            get all of the genomic features (probably genes) within a
            certain location span.

        **Arguments**
            loc
                either a location() or a cooercable location().

        **Returns**
            A vanilla list containing a bunch of keys describing any features
            sitting in the genomic span specified by location.
        """
        assert loc, "no location provided"

        try:
            loc = location(loc=loc)
        except:
            raise AssertionError, "cannot cooerce location into correct form. Location is mangled?"

        ret = []
        if loc["chr"] in self.dataByChr:
            for item in self.dataByChr[loc["chr"]]:
                #print location["left"], location["right"], item["loc"]["left"], item["loc"]["right"]
                if utils.qcollide(loc["left"], loc["right"], item["loc"]["left"], item["loc"]["right"]):
                    # make a suitable draw object
                    if 'type' not in item:
                        item["type"] = "gene" # set the type flag for gDraw
                    ret.append(item)
        return(ret)
                
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
        #if not isinstance(loc, location):
        #    loc = location(loc=loc)

        #if len(loc["chr"]) < 30: # small security measure.
        table_name = "chr_%s" % loc["chr"]

        #result = self._connection.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
        #    (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))
        
        # This is the code used in location.collide():
        #self["right"] >= loc["left"] and self["left"] <= loc["right"]
        result = self._connection.execute("SELECT left, right, strand FROM %s WHERE (right >= ? AND left <= ?)" % table_name,
            (loc["left"], loc["right"]))
    
        #result = None       
        result = result.fetchall() # safer for empty lists and reusing the cursor
        
        if result and strand: # sort out only this strand
            if strand in positive_strand_labels:
                strand_to_get = positive_strand_labels
            elif strand in negative_strand_labels:
                strand_to_get = negative_strand_labels
            
            newl = []
            for r in result:
                if r[2] in strand_to_get:
                    newl.append(r)
            result = newl                    

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
        if not self._c:
            self._c = self._connection.cursor()

        table_name = "chr_%s" % loc["chr"]

        self._c.execute("SELECT left, right, strand FROM %s WHERE (right >= ? AND left <= ?)" % table_name,
            (loc["left"], loc["right"]))
            
        return(len(self._c.fetchall()))

    def get_chromosome_names(self):
        """
        **Purpose**
            Return a list of all the valid chromosome names in the database
            
        **Arguments**
            None
            
        **Returns**
            A list of strings containing all of the chromosome names in the track
        """
        if not self._c:
            self._c = self._connection.cursor()
        
        self._c.execute("SELECT chromosome FROM main")
        r = [i[0] for i in self._c.fetchall()]
        return(set(r))

    def get_numreads_on_chromosome(self, name):
        """
        **Purpose**
            Return the number of reads on chromosme name
            
        **Arguments**
            name (Required)
                get the number of reads on chromsome 'name'
            
        **Returns**
            An integer containing the number of reads
        """
        if not self._c:
            self._c = self._connection.cursor()

        self._c.execute("SELECT chromosome, seq_reads FROM main WHERE chromosome=?", (str(name), ))
        r = self._c.fetchone()
        return(r[1])

    def get_total_num_reads(self):
        """
        **Purpose**
            Return the number total number of reads for this track.
            
        **Arguments**
            None
            
        **Returns**
            An integer containing the number of reads
        """
        if not self._c:
            self._c = self._connection.cursor()

        self._c.execute("SELECT chromosome, seq_reads FROM main")
        r = [int(i[1]) for i in self._c.fetchall()]
        
        return(sum(r))

    def _debug__print_all_tables(self):
        c = self._connection.cursor()
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
