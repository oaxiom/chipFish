"""
track, part of glbase

"""

from __future__ import division

import cPickle, sys, os, struct, math, sqlite3, zlib, time

from location import location
import utils, config
from data import positive_strand_labels, negative_strand_labels
from draw import draw
import matplotlib.cm as cm
import matplotlib.pyplot as plot
import scipy.stats as stats

from numpy import array, zeros, set_printoptions, int32, append, linspace

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
            self.meta_data = {"name": self.name,
                "source_filename": filename,
                "creation_date": time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()),
                "version": "1.01"}
            self.__setup_meta_data(self.meta_data)
        else:
            self.__load_tables(filename)
            self.__load_meta_data()

        self.__c = None
        self.__draw = None # Lazy set-up in track. Don't init draw unless needed.

    def __getitem__(self, key):
        """
        make meta data accesible like a dict
        """
        if key == "info": # catch this special key
            for k in self.meta_data:
                print "%s\t:\t%s" % (k, self.meta_data[k])
        else:
            assert key in self.meta_data, "'%s' not found in this track" % key
            return(self.meta_data[key])

    def __load_meta_data(self):
        """
        retrieve the meta data from the
        """
        c = self.__connection.cursor()

        c.execute("SELECT * FROM info")

        result = c.fetchall()

        self.meta_data = {}
        for item in result:
            self.meta_data[item[0]] = item[1]

        self.__connection.commit()
        c.close()

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
        c.execute("CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT)")

        self.__connection.commit()
        c.close()

    def __setup_meta_data(self, info_dictionary):
        """
        Load a dictionary of meta data into the info table
        """
        c = self.__connection.cursor()

        for key in info_dictionary:
            c.execute("INSERT INTO info VALUES (?, ?)", (key, info_dictionary[key]))

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

    def get_array(self, loc, strand=None, resolution=1, read_extend=0, kde_smooth=False, 
    	view_wid=0, **kargs):
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
            
            kde_smooth (Experimental)
                perform kde smoothng on the data, using the integer specified as an option.
                In this case the read_extend acts as a tag shift instead of a read_extend
                Hence set that to half of the expected shear size.

        **Returns**
            an 'numpy.array([0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
        # check if the array is already on the cache.
        # um.... 
        # Yes...
        extended_loc = location(loc=str(loc)).expand(read_extend)
        locs_required = [extended_loc] # make sure to include the read_extend, but don't modify the location to get.

        # get the reads only in the ranges required.
        # at the moment this is not implemented.
        reads = []
        for l in locs_required:
            reads += self.get_reads(l, strand)

        if kde_smooth:
            return(self.__kde_smooth(loc, reads, resolution, 0, view_wid, read_extend))

        # make a single array
        a = zeros(int( (loc["right"]-loc["left"]+resolution)/resolution ), dtype=int32)

        # work out and standardise strand option so that all forms of strand are recognised
        if strand:
            if strand in positive_strand_labels:
                strand = positive_strand_labels
            elif strand in negative_strand_labels:
                strand = negative_strand_labels

        # Bound the resolution to 1bp if int(resolution) reports <1
        min_res = int(resolution)
        if min_res == 0:
            min_res = 1

        for r in reads:
            if not strand or strand and r[2] in strand:
                if r[2] in positive_strand_labels:
                    read_left = r[0]
                    read_right = r[1] + read_extend + 1 # coords are open
                elif r[2] in negative_strand_labels:
                    read_left = r[0] - read_extend
                    read_right = r[1] + 1 # coords are open 

                rel_array_left = int((read_left - loc["left"]) / resolution)
                rel_array_right = int((read_right - loc["left"]) / resolution)

                for array_relative_location in xrange(rel_array_left, rel_array_right, 1):
                    if array_relative_location >= 0 and array_relative_location < len(a): # within array
                        a[array_relative_location] += 1

        return(a)

    def __kde_smooth(self, loc, reads, resolution, bandwidth, view_wid, read_shift=100):
        """
        Internal abstraction for kde smoothing
        
        Returns a new array 
        """
        # Shift the reads
        newr = []
        for r in reads:
            if r[2] in positive_strand_labels:
                newr.append(float((r[0] + read_shift) - loc["left"])) # Must be floats for kde to work
            elif r[2] in negative_strand_labels:
                newr.append(float((r[1] - read_shift) - loc["left"]))
        
        a = linspace(0, loc["right"] - loc["left"], view_wid)
        
        # Hack gaussian_kde()
        def covariance_factor(self):
            return 0.01
        
        kde = stats.gaussian_kde(newr)
        setattr(kde, 'covariance_factor', covariance_factor.__get__(kde, type(kde)))
        kde._compute_covariance()

        kk = kde.evaluate(a) * 1000000 # resacle to get in integer range.
        res = array(kk, dtype=int32)
        
        return(res)
                

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
            an 'numpy.array([0, 1, 2 ... n], dtype=int32)' contiginous array of integers
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
        a = zeros(int(right_most+int(resolution), int(resolution)), dtype=int32)

        for r in reads:
            # read_extend
            if r[2] in positive_strand_labels:
                left = r[0]
                right = r[1] + read_extend + 1
            elif r[2] in negative_strand_labels:
                left = r[0] - read_extend
                right = r[1] + 1

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

    def draw_pileup(self, genelist=None, key="loc", filename=None, heatmap_filename=None, bin_size=500, window_size=5000, read_extend=200, use_kde=False, **kargs): 
        """
        **Purpose**
            draw cumulative 'pileups' of the tag density in the track based on a genelist
            containing a "loc" tag
        
        **Arguments**
            genelist (Required)
                A genelist-like object containing a "loc"-like key tag
            
            key (Optional, default="key")
                the loc key to use for the locations.
                
            filename (Required)
                the name of the image file to save the pileup graph to.
                
            heatmap_filename (Optional)
                draw a heatmap, where each row is a single gene. And red is the tag density
                
            bin_size (Optional, default=500)
                bin_size to use for the heatmap
                
            window_size (Optional, default=5000)
                The window size +- around the centre of the peak to build the tag density map
                from
                
            read_extend (Optional, default=200)
                extend the read x bp either 5' or 3' depending upon the strand of the read.
                
            use_kde 
            	Use KDE versions of the tracks instead
                
        **Returns**
            If succesful returns a list containing the amalgamated pileup.
        """
        assert genelist, "you must provide a valid genelist-like object for draw_pileups()"
        assert key in genelist.linearData[0], "the genelist appears to be lacking a loc key"
        assert filename, "a filename is required for the pileup graph"
        
        n = 0
        h = 0
        pileup = None
        sampled = False
        binned_data = None
        setup_bins = False

        for i, g in enumerate(genelist):
                l = g[key].pointify().expand(window_size)
                if not use_kde:
                	a = self.get_array(l, read_extend=read_extend) # read_extend=shear size
                if use_kde:
                	a = self.get_array(l, read_extend=read_extend, kde_smooth=True, view_wid=window_size) # read_extend=tag shift
                	
                if not sampled:
                        pileup = a
                        sampled = True
                else:
                        pileup = pileup + a

                if not setup_bins:
                    binned_data = array([utils.bin_data(a, bin_size)])
                    setup_bins = True
                else:
                    #print binned_data, [utils.bin_data(a, bin_size)]
                    binned_data = append(binned_data, [utils.bin_data(a, bin_size)], axis=0)

                n += 1
                if n> 10000:
                        h += 1
                        n = 0
                        config.log.info("Done %s0,000" % h)
                        #break

        # matplotlib pileup graph
        plot.cla()
        fig = plot.figure()
        ax = fig.add_subplot(111)
        ax.plot(pileup)
        fig.savefig(filename)
        config.log.info("Saved pileup tag desity to '%s'" % filename)
        
        # Heatmap
        if not self.__draw:
            self.__draw = draw(self)
        
        if heatmap_filename:
            if "vmax" in kargs and kargs["vmax"]: 
                real_filename = self.__draw._simple_heatmap(data=binned_data, filename=heatmap_filename, dpi=150, figsize=(6, 24), vmin=0, vmax=kargs["vmax"])
            else:
                real_filename = self.__draw._simple_heatmap(data=binned_data, filename=heatmap_filename, dpi=150, figsize=(6, 24))
            config.log.info("saved heatmap to '%s'" % real_filename)

        return(pileup)
                