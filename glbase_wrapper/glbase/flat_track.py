"""
track, part of chipFish

(c) 2008-2011 oAxiom

Not for distribution.

"""

from __future__ import division

import cPickle, sys, os, struct, ConfigParser, math, sqlite3, zlib # renamed to configparser in >2.6

from location import location
from data import positive_strand_labels, negative_strand_labels

from array import array
from base_track import base_track
from draw import draw
from progress import progressbar
import config

import numpy
import matplotlib.pyplot as plot

TRACK_BLOCK_SIZE = 2000 # should go in opt, later
CACHE_SIZE = 100000 # maximum number of blocks to keep in memory.
# these are equivalent to about 800 Mb's of memory on a Fedora box
# On a 2Gb machine it is still usable, less may be a problem.
# needs to be tested on a Windows and MACOSX machine.

class flat_track(base_track):
    def __init__(self, name=None, new=False, filename=None, bin_format=None):
        """
        **Purpose**
            track definition, used for things like sequence reads across the genome
        
        **Arguments**
            name (string)
                name for the track, if not specified, it takes it from the meta_data
                If a new bd then it defaults to filename
    
            filename (string)
                directory location of the track file.
                only respected if dir_name is set.
                
            bin_format (Optional, default="i")
                the format to store the data in. 
                This is the same format as the python array, so "i" = integer,
                "f" = float
    
        """
        base_track.__init__(self, name, new, filename)
        
        if bin_format: # Change the bin_format of the db.
            self.meta_data["bin_format"] = bin_format
            self._save_meta_data()
        else:
            self.bin_format = self.meta_data["bin_format"] # I unload this one 
            
        self.bin_format = bin_format
        self.bin_len = struct.calcsize(self.bin_format)
        self.block_size = TRACK_BLOCK_SIZE # size of the blocks
        
        if new:
            self.__setup_tables(filename)
        
        # internals
        self.cache = {}
        self.cacheQ = []
        
        # Update any metadata etc, primarily the bin_len
        
    def __repr__(self):
        return("glbase.flat_track")

    def __setup_tables(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # If I make it here then base_genelist has made self._connection valid.

        c = self._connection.cursor()

        c.execute("CREATE TABLE data (blockID TEXT PRIMARY KEY, array TEXT)")

        self._connection.commit()
        c.close()

    def add_read(self, loc, strand="+", increment=1):
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

        blocks_required = ["%s:%s" % (loc["chr"], b) for b in xrange(left_most_block * self.block_size, right_most_block * self.block_size, self.block_size)]

        for blockID in blocks_required:
            # this is the span location of the block
            #block_loc = location(chr=blockID.split(":")[0], left=blockID.split(":")[1], right=int(blockID.split(":")[1])+self.block_size-1)

            # check if the block already exists
            if not self.__has_block(blockID): # not present, make a new one.
                self.__new_block(blockID)
            else:
                if not blockID in self.cacheQ: # not on cache, add it;
                    self.cacheQ.insert(0, blockID) # put the ID at the front.
                    self.cache[blockID] = self.__get_block(blockID)
            # block should now be on cache and accesable.

            bleft = int(blockID.split(":")[1])
            lleft = int(loc["left"])
            lright = int(loc["right"])
            # modify the data
            for pos in xrange(self.block_size): # iterate through the array.
                local_pos = bleft + pos # only require "left"
                # some funny stuff here:
                # right is inc'd by 1
                # as that is what is usually expected from 10-20.
                # (i.e. coords are NOT 0-based and are closed).
                if local_pos >= lleft and local_pos <= lright: # within the span to increment.
                    self.cache[blockID][pos] += increment

            self.__flush_cache()
        return(True)

    def add_score(self, loc=None, chromosome=None, left=None, right=None, score=0, **kargs):
        """
        **Purpose**
            adds a already known score at a particular location (or span). This will overwrite any previously
            exisiting value stored in the db.
            
        **Arguments**
            loc
                location.
            
            score 
                the value to insert at that location
                
        **Returns**
            Nothing
        """
        if loc:
            chrom = loc["chr"]
            left = loc["left"]
            right = loc["right"]
        else:
            chrom = chromosome    
        lleft = int(left)
        lright = int(right)
        
        left_most_block = int(abs(math.floor(left / self.block_size)))
        right_most_block = int(abs(math.ceil((right+1) / self.block_size)))

        blocks_required = ["%s:%s" % (chrom, b) for b in xrange(left_most_block * self.block_size, right_most_block * self.block_size, self.block_size)]

        for blockID in blocks_required:
            # this is the span location of the block
            #block_loc = location(chr=blockID.split(":")[0], left=blockID.split(":")[1], right=int(blockID.split(":")[1])+self.block_size-1)

            if blockID not in self.cacheQ:  
                if not self.__has_block(blockID): # A db hit, but only a check
                    self.__new_block(blockID) # no db hit
                else: # retrieve from db, big db hit
                    if not blockID in self.cacheQ: # not on cache, get it;
                        self.cacheQ.insert(0, blockID) # put the ID at the front.
                        self.cache[blockID] = self.__get_block(blockID)
            
            # block should now be on cache and accesable.
            bleft = int(blockID.split(":")[1])
            bright = bleft + TRACK_BLOCK_SIZE
            # modify the data
            #print blockID, lleft, lright, bleft, bright, score,
            for pos in xrange(lleft, lright): # iterate through the array.
                local_pos = pos - bleft # only require "left"
                # some funny stuff here:
                # right is inc'd by 1
                # as that is what is usually expected from 10-20.
                # (i.e. coords are NOT 0-based and are closed).
                if local_pos >= 0 and local_pos < TRACK_BLOCK_SIZE: # within the span to increment.
                    self.cache[blockID][local_pos] = score
                    #print local_pos, local_pos >= 0 and local_pos <= TRACK_BLOCK_SIZE
            #print self.cache[blockID]
        self.__flush_cache() # I can go over the CACHE in this routine.
        # But putting this here means I don't have to hit the db every blockID
        # Should help speed where I use a lot of new blocks.
        return(None)

    def __has_block(self, blockID):
        """
        checks if the data base has that block already.
        returns only True or False,
        does not return the block, you must use __get_block()
        to get the actual block.
        """
        if self.cache.has_key(blockID):
            return(True) # on cache, so must exist.

        c = self._connection.cursor()
        c.execute("SELECT blockID FROM data WHERE blockID=?", (blockID, ))
        result = c.fetchone()
        c.close()
        if result:
            return(True)
        return(False)

    def __commit_block(self, blockID, data):
        """
        update the block with new data.
        commit block to db.
        """
        # see if tyhe block is in the database:
        c = self._connection.cursor()
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
            data = array(self.bin_format, [0 for x in xrange(self.block_size)]) # Numpy arrays may be faster here.

        if not blockID in self.cacheQ: # not on cache
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
        if self.cache.has_key(blockID):
            return(self.cache[blockID])

        # not on the cache. get the block and put it on the cache.
        c = self._connection.cursor()
        c.execute("SELECT array FROM data WHERE blockID=?", (blockID, ))
        result = c.fetchone()
        c.close()

        if result:
            return(self.__unformat_data(result[0]))
        else:
            raise Exception, "No Block!"

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
        a = array(self.bin_format)
        a.fromstring(zlib.decompress(data))
        return(a)

    def get(self, loc, strand="+", mask_zero=False, **kargs):
        """
        **Purpose**
            get the data between location 'loc'

        **Arguments**
            loc (Required)
                a valid location or string location.

            strand (Optional, default = "+")
                strand, but only valid for stranded tracks
                
            mask_zero (Optional, defaule=False)
                return a masked numpy array with zeros masked out.

        **Returns**
            an 'array('i', [0, 1, 2 ... n])' contiginous array
        """
        try:
            if loc["chr"]: pass
        except TypeError: # probably a location string. try to cooerce
            loc = location(loc=loc)
            # Don't catch any exceptions here. Should break.

        left_most_block = int(abs(math.floor(loc["left"] / self.block_size)))
        right_most_block = int(abs(math.ceil((loc["right"]+1) / self.block_size)))

        blocks_required = ["%s:%s" % (loc["chr"], b) for b in xrange(left_most_block * self.block_size, right_most_block * self.block_size, self.block_size)]

        ret_array = array(self.bin_format, [])

        for blockID in blocks_required:
            # this is the span location of the block
            block_loc = location(chr=blockID.split(":")[0], left=blockID.split(":")[1], right=int(blockID.split(":")[1])+self.block_size-1)

            # check if the block already exists
            if self.__has_block(blockID): # it does, get it.
                block = self.__get_block(blockID)
                this_block_array_data = block # get back the usable data
            else: # block not in db, fake a block instead.
                this_block_array_data = array(self.bin_format, [0 for x in xrange(self.block_size)])

            #print "b", block
            #print self.__get_block(blockID)

            # modify the data
            for pos in xrange(self.block_size): # iterate through the array.
                local_pos = block_loc["left"] + pos
                # see add_location for details
                if local_pos < (loc["right"]+1): # still within block
                    if local_pos >= loc["left"] and local_pos <= (loc["right"]+1): # within the span to increment.
                        ret_array.append(this_block_array_data[pos])
        if mask_zero:
            mask = []
            for dd in ret_array:
                if int(dd*10000000) == 0: # delete if 10 sig figs close to 0
                    mask.append(1)
                else:
                    mask.append(0)
            # This may produce a warning, but apparently I can safely ignore it
            ret_array = numpy.ma.masked_array(ret_array, mask=mask)
        return(ret_array)

    def finalise(self):
        """
        {Override)
        I have to override the base_track class
        finalise the database (shrink unused edit space)
        dump useless bits etc.
        You must call this! to finalise the db.
        Or get() will not work!
        This copies the cache onto disk and closes the db.
        """
        for blockID in self.cache:
            self.__commit_block(blockID, self.cache[blockID])
        self.cache = {} # Kill caches.
        self.cacheQ = []
        self._connection.commit() # commit all the __commit_block()

        base_track.finalise(self)
        
    def pileup(self, genelists=None, filename=None, bandwidth=None, average=True, 
        background=None, mask_zero=False, respect_strand=True, **kargs):
        """
        **Purpose**
            Draw a set of pileup count scores (averages or cumulative scores)
        
        **Arguments**
            genelists
                A genelist with a "loc" key
            
            filename
                The filename to save the image to
                
            bandwidth (Optional, default=None)
                the number of base pairs to use around the centre of the location
                If set to None then it will use the location as specified. 
        
            average (Optional, default=True)
                use the average score if set to True (i.e. divide by the number of items)
                Or use the cumulative score (the total) if set to False
                
            background (Optional)
                You can supply a list of background coordinates if you want. The lists
                must contain a "loc" key.
                
            mask_zero (Optional, default=False)
                flat_tracks are continuous and must have a value at all base pairs. I fill
                that value with 0
                However, in some flat_tracks 0 actually means 'no data' and I want to ignore 
                those values in the pileup. If that is the case for you, set mask_zero to
                True.
                
            <other graph-related args>
            pileup also respects:
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title

            respect_strand (Optional, default=False)
                If available, respect the orientation of the strand from the genelist.
                This is useful if you are, say, using the TSS's and want to maintain the
                orientation with respect to the transcription direction.
                
        **Returns**
            (data, back)
            The data from the line graph. 
            back will be the average of the list of background peaks, but data
            will be the last entry in the peaklists (if a list of peaks) or will correspond to the
            only peaklist provided. e.g.:
            
            data, back = pileup(genelists=[data1, data2, THISDATAWILLBERETURNED] ...)
            
            or:
            
            data, back = pileup(genelists=THISDATAWILLBERETURNED ...)
        """
        assert genelists, "genelists is None?"
        
        if not isinstance(genelists, list):
            genelists = [genelists] # make a one item'd list
    
        if background:
            if not isinstance(background, list):
                background = [background] # make a one item'd list
    
        # flats have lazy setup of draw:
        if not self._draw:
            self._draw = draw(self)
        
        fig = self._draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        x = None
        if bandwidth:
            x = numpy.arange(bandwidth*2) - int((bandwidth*2)/2)
        
        if bandwidth:
            loc_span = bandwidth*2
        else:
            loc_span = len(genelists[0].linearData[0]["loc"]) # I have to assume all locs are identical.
        
        for gl in genelists:
            if bandwidth:
                hist = numpy.zeros(bandwidth*2)
                counts = numpy.zeros(bandwidth*2)
            else:
                x = numpy.arange(loc_span) - int(loc_span/2)
                hist = numpy.zeros(loc_span)
                counts = numpy.zeros(loc_span) # used to get the average.
                
            for i in gl:
                if bandwidth:
                    l = i["loc"].pointify()
                    l = l.expand(bandwidth)
                    a = self.get(l)[0:bandwidth*2] # mask_zero is NOT asked of here. because I need to ignore zeros for the average calculation (below)
                else:
                    a = self.get(i["loc"])[0:loc_span]
                
                if respect_strand:
                    # positive strand is always correct, so I leave as is.
                    # For the reverse strand all I have to do is flip the array.
                    if i["strand"] in negative_strand_labels:
                        a = a[::-1]
                
                if a: # It's possible that get() will return nothing 
                    # For example if you send bad chromosome names or the locations are nonsensical (things like:
                    # chr9_GL000201_RANDOM:-500-1500
                    hist += a

                if mask_zero: # surely a better way of doing this...
                    t = numpy.zeros(loc_span)
                    for ee, xx in enumerate(a):
                        if xx > 0:
                            t[ee] = 1.0
                    counts += t                 

            if average and mask_zero:
                hist /= counts
            elif average and not mask_zero:
                hist /= len(gl)
            
            ax.plot(x, hist, label=gl.name, alpha=0.7)
            
        bkgd = None
        if background:
            if bandwidth:
                bkgd = numpy.zeros(bandwidth*2)
                counts = numpy.zeros(bandwidth*2)
            else:
                bkgd = numpy.zeros(loc_span)
                counts = numpy.zeros(loc_span)
                
            bkgd_items = 0
            p = progressbar(len(background))
            for i, back in enumerate(background):
                for b in back:
                    if bandwidth:
                        l = b["loc"].pointify()
                        l = l.expand(bandwidth)
                        a = self.get(l)[0:bandwidth*2]
                    else:
                        a = self.get(b["loc"])[0:loc_span]
                    bkgd_items += 1
                    
                    if respect_strand:
                        # positive strand is always correct, so I leave as is.
                        # For the reverse strand all I have to do is flip the array.
                        if b["strand"] in negative_strand_labels:
                            a = a[::-1]
                
                    bkgd += a
                    if mask_zero:
                        t = numpy.zeros(loc_span)
                        for ee, xx in enumerate(a):
                            if xx > 0:
                                t[ee] = 1.0
                        counts += t                 

                if average and mask_zero:
                    bkgd /= counts
                elif average and not mask_zero:
                    bkgd /= bkgd_items
                
                if i == 0: # add only a single legend.
                    ax.plot(x, bkgd, color="grey", alpha=0.3, label="Random Background") 
                else:
                    ax.plot(x, bkgd, color="grey", alpha=0.3) 
                
                # reset arrays
                bkgd = numpy.zeros(len(bkgd))
                counts = numpy.zeros(len(counts))
            
                p.update(i)

        else:
            bkgd = None
        
        leg = ax.legend()
        [t.set_fontsize(7) for t in leg.get_texts()]
        ax.set_ylabel("Magnitude")
        ax.set_xlabel("Base pairs around centre (bp)")
        
        self._draw.do_common_args(ax, **kargs)
               
        actual_filename = self._draw.savefigure(fig, filename)
        
        config.log.info("pileup(): Saved '%s'" % actual_filename)
        return(hist, bkgd)