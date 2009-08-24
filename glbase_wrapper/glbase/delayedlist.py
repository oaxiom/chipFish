"""
Think of me as a delayed version of geneList

* the returned arrangement of
* this is the only special case of geneList (delayed)
* should I make this part of genelist?
* I think this can be generalised instead to be a 'delayed' genelist, with
  the same interface as genelist (where possible) but the data remains on the disk.
  See the Klf_Snipper for an example where I use this

"""
import sys, os, time, copy

import config
import utils
from flags import *
from array import array as qarray
from draw import draw
from taglist import taglist
from genelist import genelist
from history import historyContainer
from errors import AssertionError

# try to load non-standard libs.
try:
    import matplotlib.pyplot as plot
    MATPLOTLIB_AVAIL = True
except:
    print "Warning: matplotlib not available or not installed"
    MATPLOTLIB_AVAIL = False

class delayedlist(genelist):
    """
    delayedlist is a descendent of genelist - except the data never
    makes it into memory.
    """
    def __init__(self, filename=None, **kargs):
        genelist.__init__(self) # no kargs though.

        assert filename, "No Filename"
        assert os.path.exists(filename), "%s not found" % (filename)
        assert kargs.has_key("format"), "You must provide a format for delayedlist. I cannot guess its format."

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullpath = filename
        self.filehandle = None

        self.format = kargs["format"] # override default
        self._optimiseData()

    def collide(self, **kargs):
        """
        Note: the 'logic' command is not supported for delayedlists
        only "and" can be performed.
        """
        self._optimiseData()

        assert kargs.has_key("peaklist") or kargs.has_key("genelist") or kargs.has_key("microarray"), "You must provide a genelist-like object"
        assert kargs.has_key("loc_key"), "You must provide a 'loc_key' name"
        if kargs.has_key("peaklist"): gene_list = kargs["peaklist"]
        elif kargs.has_key("genelist"): gene_list = kargs["genelist"]
        elif kargs.has_key("microarray"): gene_list = kargs["microarray"]
        assert gene_list[0].has_key(kargs["loc_key"])
        assert self.__iter__().next().has_key(kargs["loc_key"]) # get an item and test it
        self._optimiseData()

        delta = 200
        if kargs.has_key("delta"): delta = kargs["delta"]

        return(genelist.collide(self, genelist=gene_list, loc_key=kargs["loc_key"], delta=delta, merge=True))

    def overlap(self, **kargs):
        """
        Note: the 'logic' command is not supported for delayedlists
        only "and" can be performed.
        """
        self._optimiseData()

        assert kargs.has_key("peaklist") or kargs.has_key("genelist") or kargs.has_key("microarray"), "You must provide a genelist-like object"
        assert kargs.has_key("loc_key"), "You must provide a 'loc_key' name"
        if kargs.has_key("peaklist"): gene_list = kargs["peaklist"]
        elif kargs.has_key("genelist"): gene_list = kargs["genelist"]
        elif kargs.has_key("microarray"): gene_list = kargs["microarray"]
        assert gene_list[0].has_key(kargs["loc_key"])
        assert self.__iter__().next().has_key(kargs["loc_key"]) # get an item and test it
        self._optimiseData()

        delta = 200
        if kargs.has_key("delta"): delta = kargs["delta"]

        return(genelist.overlap(self, genelist=gene_list, loc_key=kargs["loc_key"], delta=delta, merge=True))

    def __len__(self):
        return(-1) # send back a dummy length, even though it's meaningless.

    def __getitem__(self, index):
        """
        (Override)
        confers a = geneList[0] behaviour
        This is mostly broken. It returns only the very first entry,
        whatever 'index' is sent to it, it disregards it.
        This only continues to exist for compatability with some internal
        routines.
        """
        self._optimiseData()
        return(self.__iter__().next())
        self._optimiseData()

    def __iter__(self):
        """
        (Override)
        make the geneList behave like a normal iterator (list)
        """
        try:
            for index, column in enumerate(self.linearData):
                # format the data to look like a genuine entry.
                #if self.format.has_key("debug") and self.format["debug"]:
                #    print "%s:'%s'" % (index, column)
                #    if isinstance(self.format["debug"], int) and index > self.format["debug"]: break # If an integer, collect that many items.

                # skiplines is implemented in _optimiseData() (i.e. open and reset the file).

                d = {}
                if column: # list is empty, so omit.
                    if (not (column[0] in typical_headers)):
                        for key in self.format:
                            if not (key in ignorekeys): # ignore these tags
                                if not d.has_key(key):
                                    d[key] = {}
                                if isinstance(self.format[key], dict) and (self.format[key].has_key("code")):
                                    # a code block insertion goes here - any valid lib and one line python code fragment
                                    # store it as a dict with the key "code"
                                    d[key] = eval(self.format[key]["code"])
                                else:
                                    d[key] = self._guessDataType(column[self.format[key]])
                yield d
        except:
            raise StopIteration

    def _optimiseData(self):
        """
        (Override)
        Impossible to optimise the data.
        so we just reset the entry point.
        """
        if self.filehandle: self.filehandle.close()

        self.filehandle = open(self.fullpath, "rU")
        if self.format.has_key("dialect"):
            self.linearData = csv.reader(self.filehandle, dialect=self.format["dialect"])
        else:
            self.linearData = csv.reader(self.filehandle)

        if self.format.has_key("skiplines"):
            if self.format["skiplines"] != -1: # no skuipped klines, good to go.
                for i, x in enumerate(self.linearData):
                    if i == self.format["skiplines"]:
                        break
        return(True)

    def __str__(self):
        self._optimiseData()

        # load in a bunch of data and dump it into self.linearData
        temp_data = []
        for index, item in enumerate(self):
            temp_data.append(item)
            if index > config.NUM_ITEMS_TO_PRINT-2:
                break # only get the first n data.s
        self.linearData = temp_data
        ret = genelist.__str__(self)
        self._optimiseData()
        return("%s\r\nOnly the first %s entries are shown" %(ret, config.NUM_ITEMS_TO_PRINT-1))

    def save(self):
        print "Error: Cannot save a binary representation of a delayedlist"

    def saveCSV(self, **kargs):
        print "Error: delayedlists do not support saveCSV()"

    def getChIPSeqTags(self, gene_list, bins=[x for x in xrange(-5000, 5000, 100)], bSaveMergedImages=True):
        """
        **Purpose**

        Count the number of chip seq tags that lie under some arrangement of the 'bins'

        **Arguments**

        gene_list
            your gene list of variant

        bins
            an array of bins spannign the range of base pairs you want to bin.
            e.g. [x for x in range(-5000, 5000, 100)]

        bSaveMergedImages
            undocumented.

        **Result**

        returns a taglist list (a descendent of geneList)
        that behaves very much like gene list.
        However, any sepcial methods from the input gene_list will be lost
        in preference to taglist's methods.
        the chip_seq tag frequency is stored in the "chip_tags" key
        """
        print "Warning: getChIPSeqTags() is slow (depending upon the size of your tag file)"
        print "         press Ctrl+C to interupt."

        # test gene_list has a valid sorted
        try:
            gene_list.dataByChr["1"]
        except:
            print "Error: List does not have a valid location tag"
            return(False)

        try:
            gene_list.dataByChr["1"][0]["tss_loc"]
        except:
            print "Error: List does not have a valid tss_loc tag"
            return(False)

        newl = taglist(bins, gene_list.name) # loses any non-standard methods... :(
        newl.linearData = []

        #def mapAgainstSeqFile(path, seq_filename, listofgenes, bins=[-5000, -4000, -3000, -2000, -1000, 0, 1000, 2000, 3000, 4000, 5000], locCol=0, symmetric=False):
        oh = open(self.fullpath, "rU")

        results = [] # fill a blank array to store the scores.
        for item in xrange(len(gene_list)+1):
            results.append(qarray("L", [0 for item in bins]))
        flatBin = qarray("L", [0 for item in bins])

        binLeft = bins[0]
        binRight = bins[-1]
        n = 0 # counters
        t = 0
        f = 0

        start = time.time()

        for line in oh:
            s = line.split("\t") # avoid the overhead of csv.reader.
            n += 1
            chr = s[self.format["chr"]].replace("chr", "")
            left = int(s[self.format["left"]])
            right = int(s[self.format["right"]])
            if gene_list.isChromosomeAvailable(chr):
                for data in gene_list.dataByChr[chr]:
                    loc = utils.getLocation(data["tss_loc"])
                    #print data

                    tss = loc["left"] # always for tss_loc

                    mid_tag = left + ((right - left)/2)
                    #print tss, binLeft, mid_tag, binRight
                    if ((tss+binLeft) < mid_tag) and ((tss+binRight) > mid_tag):
                        if data["strand"] == "+":
                            local_tag_loc = (mid_tag + binLeft) - (tss + binLeft)
                        elif data["strand"] == "-":
                            local_tag_loc = (tss + binLeft) - (mid_tag + binLeft)
                        #print local_tag_loc, mid_tag

                        l = bins[0]
                        for i, b in enumerate(bins[1:]):
                            #print "l:%s b:%s" % (l,b)
                            if (local_tag_loc > l and local_tag_loc < b):
                                binSet = results[data["n"]] # get the location in the original list;
                                binSet[i] += 1
                                flatBin[i] += 1
                                f += 1
                                #print "g:", tss, item[1], "t:", mid_tag, "b:", bins[i], "-", bins[i+1]
                                #break
                            l = b

            if n > 1000000: # counter
                t += 1
                print "Info: Done: %s,000,000 - Found: %s tags" % (t, f)
                n = 0
                #break

        oh.close() # finished with the file.

        # load my newl
        for index, item in enumerate(gene_list.linearData):
            c = copy.deepcopy(item) # if you don't copy, it will mangle the original list...
            if not c.has_key("chip_tags"):
                c["chip_tags"] = {}
            c["chip_tags"][self.name] = results[index]
            newl.linearData.append(c)

        # draw some pictures of flatBin.

        end = time.time()
        print "Info: Time taken: %i mins" % ((end - start) / 60)

        print "Info: All done, found:", f
        if MATPLOTLIB_AVAIL:
            self._drawMerged(flatBin, gene_list.path)

        newl._optimiseData()
        return(newl)

    def _drawMerged(self, flatBin, path, window=1):
        """
        (Internal)
        """
        plot.cla()
        if window > 1:
            n = utils.movingAverage(flatBin, window)
        else:
            n = flatBin
        plot.plot(n)
        plot.savefig(os.path.join(path, "ChIP_merge_%s.png" % self.name))
        print("Info: Saved a merge of the ChIP-peak to: %s" % os.path.join(path, "ChIP_merge_%s.png" % self.name))
