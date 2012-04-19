"""
track, part of glbase

"""

from __future__ import division

import cPickle, sys, os, struct, math, sqlite3, zlib, time, csv

from progress import progressbar
from errors import AssertionError
from location import location
import utils, config
from data import positive_strand_labels, negative_strand_labels
from draw import draw
import matplotlib.cm as cm
import matplotlib.pyplot as plot
import scipy.stats as stats
from base_track import base_track

import numpy
from numpy import array, zeros, set_printoptions, int32, append, linspace, argmax, amax, delete

TRACK_CACHE_SIZE = 10 # number of track segments to cache.

class track(base_track):
    """
    track definition, used for things like sequence reads across the genome
    """
    def __init__(self, name=None, new=False, filename=None, **kargs):
        """
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
        base_track.__init__(self, name, new, filename)
        
        if new:
            self.__setup_tables(filename)

    def __repr__(self):
        return("glbase.track")

    def __setup_tables(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # If i make it here, then a lot of grunt work already done in base_track
        c = self._connection.cursor()

        c.execute("CREATE TABLE main (chromosome TEXT PRIMARY KEY, seq_reads INT)")

        self._connection.commit()
        c.close()

    '''
    # No plans to implement this atm.
    def __get_bucket_id(self, chromosome, left, right):
        """
        return a list of table names (or bucketids) required for left and right
        """
        # work out which of the buckets is required:
        left_buck = int((loc["left"]-1-delta)/config.bucket_size)*config.bucket_size
        right_buck = int((loc["right"]+delta)/config.bucket_size)*config.bucket_size
        buckets_reqd = range(left_buck, right_buck+config.bucket_size, config.bucket_size) # make sure to get the right spanning and left spanning sites
        return(["%s_%s" % (chromsome, bid) for bid in buckets_reqd])
    '''

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
        c.execute("CREATE TABLE %s (left INT, right INT, strand TEXT)" % (table_name, ))

        # link new table to old table.
        # how do I do that?!?!?

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

        c = self._connection.cursor()

        # insert location into new array.
        table_name = "chr_%s" % str(loc["chr"])

        # get the old number of seq_reads
        c.execute("SELECT seq_reads FROM main WHERE chromosome=?", (loc["chr"], ))
        current_seq_reads = c.fetchone()[0] # always returns a tuple

        c.execute("UPDATE main SET seq_reads=? WHERE chromosome=?", (current_seq_reads+1, loc["chr"]))

        # add the location to the seq table:
        c.execute("INSERT INTO %s VALUES (?, ?, ?)" % table_name, (loc["left"], loc["right"], strand))

        c.close()

    def get(self, loc, resolution=1, read_extend=0, kde_smooth=False, 
        view_wid=0, **kargs):
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

        **Returns**
            an 'numpy.array([0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
           
        extended_loc = location(loc=str(loc)).expand(read_extend)
        locs_required = [extended_loc] # make sure to include the read_extend, but don't modify the location to get.

        # get the reads only in the ranges required.
        # at the moment this is not implemented.
        reads = []
        for l in locs_required:
            reads += self.get_reads(l)

        if kde_smooth:
            return(self.__kde_smooth(loc, reads, resolution, 0, view_wid, read_extend))

        # make a single array
        a = zeros(int( (loc["right"]-loc["left"]+resolution)/resolution ), dtype=int32)

        # Bound the resolution to 1bp if int(resolution) reports <1
        min_res = int(resolution)
        if min_res == 0:
            min_res = 1

        for r in self.c:
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
            return 0.02
        
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
        if strand: 
            raise NotImplementedError, "Eh... strand not supported yet..."

        if not self._c:
            self._c = self._connection.cursor()
            
        table_name = "chr_%s" % chromosome
        self._c.execute("SELECT * FROM %s" % table_name)
        reads = self._c.fetchall()

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
        if not isinstance(loc, location):
            loc = location(loc=loc)
            
        if not self._c:
            self._c = self._connection.cursor()

        if len(loc["chr"]) < 30:
            table_name = "chr_%s" % loc["chr"]

        # ~4.11 s 946 reads
        #self._c.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
        #    (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))
            
        # Is it faster to just get all the nearby reads and sort them out in Python?
        self._c.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right)" % table_name,
            (loc["left"], loc["left"], loc["right"], loc["right"]))
        
        #result = self._c.fetchall()

        #return(result)

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

        self._c.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
            (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))

        # this is here in the hope that in the future cursor.rowcount
        # will be correctly supported...
        # at the moment probably provides little or no benefit.

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
        return(r)

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

    def pileup(self, genelist=None, key="loc", filename=None, heatmap_filename=None, 
        bin_size=500, window_size=5000, read_extend=200, use_kde=False, simple_cleanup=False,
        only_do=False, stranded=False, respect_strand=True, raw_tag_filename=None,
        **kargs): 
        """
        **Purpose**
            draw cumulative 'pileups' of the tag density in the track based on a genelist
            containing a "loc" tag
        
        **Arguments**
            genelist (Required)
                A genelist-like object containing a "loc"-like key tag
            
            key (Optional, default="loc")
                the key to use for the locations. Must contain valid location data:
                chr1:1000-1011 etc. draw_pileup() will use the centre of the location if it is a 
                span.
                
            filename (Required)
                the name of the image file to save the pileup graph to.
                
            heatmap_filename (Optional)
                draw a heatmap, where each row is a single gene. And red is the tag density
            
            raw_tag_filename (Optional)
                Save a tsv file that contains the heatmap values for each row of the genelist.
            
            bin_size (Optional, default=500)
                bin_size to use for the heatmap
                
            window_size (Optional, default=5000)
                The window size +- around the centre of the peak to build the tag density map
                from
                
            read_extend (Optional, default=200)
                extend the read x bp either 5' or 3' depending upon the strand of the read.
                If use_kde is true then this will be the 'tag shift' value instead.
                
            use_kde (Optional)
                Use KDE versions of the tracks instead (Experimental)
                
            simple_cleanup (False)
                remove rows from the pileup that have < simple_cleanup tag counts
                
            stranded (Optional, default=False)
                build a stranded pileup, with + reads in blue and - reads in red
                
            respect_strand (Optional, default=True)
                If available, respect the orientation of the strand from the genelist.
                This is useful if you are, say, using the TSS's and want to maintain the
                orientation with respect to the transcription direction.
                            
        **Returns**
            If succesful returns a list of lists containing the a single entry for each
            entry in the original genelist (in the same order as the genelist).
        """
        assert filename, "you must specify a filename"
        assert genelist, "you must specify a genelist"
        assert key in genelist.linearData[0], "the genelist appears to be lacking a loc key"
        
        if stranded:
            return(self.__draw_pileup_stranded(genelist, filename, window_size, **kargs))
        else:
            return(self.__draw_pileup_normal(genelist, key, filename, heatmap_filename, 
                bin_size, window_size, read_extend, use_kde, simple_cleanup, 
                only_do=only_do, raw_tag_filename=raw_tag_filename, respect_strand=respect_strand, **kargs))
                
    def __draw_pileup_stranded(self, genelist=None, filename=None, bandwidth=300, **kargs):
        """
        **Purpose**
            Build a histogram density plot of the paired reads. 
            This will estimate the approximate observed shear size in your chip-seq data.
            It pairs the data together for all pairs of all reads within a specified bandwidth
            then outputs a plot of the resulting histogram.
            
        **Arguments**
            filename
                the filename to save the image(s) to.
                
            genelist (Required)
                some sort of genelistlike object with a 'loc' key containing genomic coordinates
            
            bandwidth (Optional, default=300)
                area around the centre of the peak to build the cumulative distributions.
        
        **Returns**
            None
            and some image files in <base_filename>.<draw_mode>
        """      
        # step along each chromosome. Quit if there are no reads in the window
        
        if not self._c:
            self._c = self._connection.cursor()
        
        hist = {"+": zeros(bandwidth+1), "-": zeros(bandwidth+1)}
        
        p = progressbar(len(genelist))
        for n, read in enumerate(genelist):
            loc = read["loc"].pointify().expand(bandwidth//2)
            for s in self.get_reads(loc): # get reads returns all overlapping reads. I need to trim 
                # the edges to stop negative array positions
                
                loc_left = s[0] - loc["left"]
                loc_right = s[1] - loc["left"]

                if s[2] == "+" and loc_left > 0:
                    loc_left = s[0] - loc["left"]
                    hist["+"][loc_left] += 1
                elif s[2] == "-" and loc_right < bandwidth-1:
                    loc_right = s[1] - loc["left"]
                    hist["-"][loc_right] += 1

            p.update(n)
            
        if not self._draw:
            self._draw = draw(self)
            
        # now plot:
        fig = self._draw.getfigure()
        ax = fig.add_subplot(111)
        ax.plot(hist["+"], color="blue")
        ax.plot(hist["-"], color="red")
        
        max_left = argmax(hist["+"])
        max_right = argmax(hist["-"])
        
        realfilename = self._draw._saveFigure(fig, filename)
        config.log.info("Saved shear_size_pileup '%s'" % realfilename)
        return(hist)
    
    def __draw_pileup_normal(self, genelist=None, key="loc", filename=None, heatmap_filename=None, 
        bin_size=500, window_size=5000, read_extend=200, use_kde=False, simple_cleanup=False,
        only_do=False, respect_strand=True, raw_tag_filename=None, mask_zero=False,
        **kargs): 
        """
        The normal pileup views
        """
        # See if there is a proper stransd key in there somewhere:
        if respect_strand and (not "strand" in genelist.linearData[0]):
            config.log.warning("I could not find the 'strand' key, setting respect_strand to False")
            respect_strand = False
            
        n = 0
        h = 0
        pileup = None
        sampled = False
        binned_data = None
        setup_bins = False

        p = progressbar(len(genelist))
        for i, g in enumerate(genelist):
            l = g[key].pointify().expand(window_size)
            if not use_kde:
                a = self.get(l, read_extend=read_extend) # read_extend=shear size
            if use_kde:
                a = self.get(l, read_extend=read_extend, kde_smooth=True, view_wid=window_size) # read_extend=tag shift
                
            if respect_strand:
                # positive strand is always correct, so I leave as is.
                # For the reverse strand all I have to do is flip the array.
                if g["strand"] in negative_strand_labels:
                    a = a[::-1]
            
            if sum(a) > simple_cleanup: # Only keep if tag count is > simple_cleanup
                if not sampled: # numpy __nonzero__ retardedness
                    pileup = a
                    sampled = True
                else:
                    pileup = pileup + a

                if heatmap_filename or raw_tag_filename:
                    if not setup_bins:
                        binned_data = array([utils.bin_data(a, bin_size)])
                        setup_bins = True
                    else:
                        #print binned_data, [utils.bin_data(a, bin_size)]
                        binned_data = append(binned_data, [utils.bin_data(a, bin_size)], axis=0)
                
            if only_do and n > only_do: 
                #print only_do, n
                break
            p.update(i)

        if not self._draw:
            self._draw = draw()

        # matplotlib pileup graph
        fig = self._draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.plot(pileup)
        real_filename = self._draw.savefigure(fig, filename)
        config.log.info("Saved pileup tag density to '%s'" % filename)
               
        other_args = {}
        if "vmax" in kargs:
            other_args["vmax"] = kargs["vmax"]
        if "vmin" in kargs:
            other_args["vmin"] = kargs["vmin"]

        # spin this off into a .heatmap() method?
        if heatmap_filename or raw_tag_filename:
            binned_data = numpy.delete(binned_data, numpy.s_[-1:], 1) # kill the rightmost empty col.
            
            if raw_tag_filename:
                binned_data = numpy.delete(binned_data, numpy.s_[-1:], 1)
                oh = open(raw_tag_filename, "w")
                writer = csv.writer(oh, dialect=csv.excel_tab)
                writer.writerows(binned_data)
                oh.close()
                config.log.info("saved raw_tag_file to '%s'" % raw_tag_filename)
                
            if heatmap_filename:
                if other_args: 
                    real_filename = self._draw._simple_heatmap(data=binned_data, filename=heatmap_filename, 
                        dpi=150, figsize=(6, 24), aspect="long", **other_args)
                else:
                    real_filename = self._draw._simple_heatmap(data=binned_data, filename=heatmap_filename, 
                        dpi=150, figsize=(6, 24), aspect="long")
                config.log.info("saved heatmap to '%s'" % real_filename)

        return(pileup)

    def finalise(self):
        """
        finalise the database (shrink unused edit space)
        dump useless bits etc.
        You must call this! to finalise the db.
        Or get() will not work!
        This copies the cache onto disk and closes the db.
        """
        # do a commit
        self._save_meta_data()
        self._connection.commit()
        self._connection.execute("VACUUM")
        self._connection.commit()

if __name__ == "__main__":
    import random
    from location import location
    from genelist import genelist
    
    t = track(filename="test.trk2", name="test", new=True)
    for n in xrange(0, 100000):
        l = random.randint(0, 10000000)
        t.add_location(location(chr="1", left=l, right=l+35), "+")
    t.finalise()
    
    # fake a bed
    newb = []
    for n in xrange(0, 100):
        l = random.randint(0, 10000000)
        newb.append({"loc": location(chr="1", left=l, right=l+35), "strand": "+"})
    bed = genelist()
    bed.load_list(newb)
    
    t = track(filename="test.trk2")
                
    import cProfile, pstats
    cProfile.run("t.pileup(genelist=bed, filename='test.png', bin_size=10, window_size=1000, raw_tag_filename='Sox2Oct4_P300.tsv')", "profile.pro")
    p = pstats.Stats("profile.pro")
    p.strip_dirs().sort_stats("time").print_stats()