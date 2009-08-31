"""
A special ordered list of genes.
behaves like a normal list, but each element contains a heterogenous set of data.

**to-do list**

* add an option for extra_header [] to format. [skiplines?]
* store more stuff in optimiseData for faster searching/indexing/matching. (pre-sorted list pointers?)
* saveCSV does not save the labels in any order, they could come out all over the place! [Fixed?]
* map is evil slow.
* get rid of this idea of 'paths' and just use filenames (I can store paths internally) [Partially done]?
* nicely clean up from errors. [Partially Done]
* logging [Started]
* uniformly support location class.
* set only copies the data - does not produce a genelist. [Think this is fixed]
* make all methods keyword args based. [Done? Not checked]
* Always inherit from the right side????
* ALL non-internal commands should use **kargs. (Only partially implemented)
* make add, sub, and, or work on non-equivalent lists. [Partially implemented]
* wider use of assertion to declare required arguments. [Completed? Implemnted?]
"""

import sys, os, csv, copy, random, cPickle

from operator import itemgetter

import config
from flags import *
from helpers import *
from location import location
from draw import draw
from history import historyContainer
from errors import AssertionError, UnRecognisedCSVFormatError, UnrecognisedFileFormatError

class genelist:
    """
    This is essentially a class container for any arrangement of heterogenous data.

    it is good for dealing with csv files with arbitrary columns - arranging
    them by keys and then allowing cross-matching and data extraction.
    geneList is the generic implementation, and you should probably avoid it as it is a bit raw
    instead use one of the wrappers, which are geneLists with more sensible defaults
    and drawing classes for known arrangements of data. Always use them in preference
    to a vanilla genelist.
    """
    def __init__(self, **kargs):
        self.linearData = []
        self.dataByChr = None # this is private, use get by loc.
        self.debug = False
        self.draw = draw(self)
        self.name = "Generic List"
        self._history = historyContainer(self) # a list of historyItem's

        format = sniffer
        if kargs.has_key("name"): self.name = kargs["name"]
        if kargs.has_key("format"): format = kargs["format"] # I expect a filename= is coming.
        if kargs.has_key("filename"):
            if kargs.has_key("format"):
                self.load(filename=kargs["filename"], format=kargs["format"])
            else:
                self.load(filename=kargs["filename"])

    def load(self, filename=None, format=None, **kargs):
        """
        **Purpose**

        load a file into the genelist. load will attempt to load the file
        based on the filename, unless a format is specified.

        **Arguments**

        filename
            absolute filename (including path) to the actual file.
            can include path short cuts (e.g. "./", "../" etc)

        format (Optional, default = "sniffer" (ie. guess))
            format specifer, see flags.py and helpers.py and the
            documentation on how to write a valid format specifier

        **Result**

        fills the genelist with the data from the file as specified by
        the format specifier.
        """
        assert filename, "No filename specified"
        assert os.path.exists(os.path.realpath(filename)), "File %s not found" % filename

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullfilename = filename
        self.name = self.filename.split(".")[0]

        if format:
            if format.has_key("special"): # special loads
                if format["special"] == "fasta":
                    self.linearData = utils.convertFASTAtoDict(filename)
                    self._optimiseData()
                    return(True)
        else:
            csv_headers = frozenset(["csv", "xls", "tsv", "txt"])
            if filename.split(".")[-1] in csv_headers: # check the last one for a csv-like header
                if not format:
                    self.loadCSV(filename=filename, format=sniffer)
                else:
                    self.loadCSV(filename=filename, format=format)
            elif filename.split(".")[-1] in ["glb"]:
                self = load(self.filename) # will this work?
            else:
                if not format:
                    self.loadCSV(filename=filename, format=sniffer)
                else:
                    self.loadCSV(filename=filename, format=format)
            return(True) # must have made it to one - if it fails it should trigger
            # UnrecognisedFileFormatError...
        # Oh dear, no recognised format.
        raise UnrecognisedFileFormatError, ("File format is not recognised", filename, format)

    def loadCSV(self, filename=None, format=sniffer, **kargs):
        """
        **Purpose**

        load a CSV file into the genelist

        **Arguments**

        filename
            absolute filename (including path) to the actual file.
            can include path short cuts (e.g. "./", "../" etc)

        format (Optional, default = "sniffer" (ie. guess))
            format specifer, see flags.py and helpers.py and the
            documentation on how to write a valid format specifier

        name (Optional, Default = based on the filename)
            name of the genelist (stored in self.name, not really used
            widely to be honest)

        **Result**

        fills the genelist with the CSV table.
        """
        assert os.path.exists(os.path.realpath(filename)), "File %s not found" % filename

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullfilename = filename
        self.name = self.filename.split(".")[0]

        # make a guess if it is a tsv vs csv

        self._loadCSV(filename=self.fullfilename, format=format, **kargs)
        try:
            pass
        except:
            try: # try again, guessing it might be a tsv
                format["dialect"] = csv.excel_tab
                self._loadCSV(filename=self.fullfilename, format=format, **kargs)
            except: # oh dear. Die.
                raise UnRecognisedCSVFormatError, ("csv format mangled, the csv does not fit the format specifier", self.fullfilename, format)

    def _loadCSV(self, **kargs):
        """
        (Internal)

        Wrapper around the actual loadCSV()
        """
        assert kargs.has_key("filename"), "No filename specified"
        assert kargs.has_key("format"), "_loadCSV requres a format specifier"

        filename = kargs["filename"]
        format = kargs["format"]

        temp_data = []
        oh = open(filename, "rU")
        if format.has_key("dialect"):
            reader = csv.reader(oh, dialect=format["dialect"])
        else:
            reader = csv.reader(oh)

        if format.has_key("sniffer"):
            # I need to construct my own format
            format = {}
            for top_line in reader:
                for index, key in enumerate(top_line): # get all the potential keys.
                    format[key] = index
                break

        if format.has_key("skiplines"):
            skiplines = format["skiplines"]
        else:
            skiplines = 0 # skip any header row by default.

        for index, column in enumerate(reader): # This is cryptically called column, when it is actually row.
            if format.has_key("debug") and format["debug"]:
                print "%s:'%s'" % (index, column)
                if isinstance(format["debug"], int) and index > format["debug"]: break # If an integer, collect that many items.

            if index > skiplines:
                # there is a reason for that, it is so that in the formats it appears:
                # "refseq": column[1] # see :)
                if column: # list is empty, so omit.
                    if (not (column[0] in typical_headers)):
                        d = {}
                        for key in format:
                            if not (key in ignorekeys): # ignore these tags
                                if not d.has_key(key):
                                    d[key] = {}
                                if isinstance(format[key], dict) and (format[key].has_key("code")):
                                    # a code block insertion goes here - any valid lib and one line python code fragment
                                    # store it as a dict with the key "code"
                                    d[key] = eval(format[key]["code"])
                                else:
                                    d[key] = self._guessDataType(column[format[key]])
                        temp_data.append(d)
        oh.close()

        d = []
        if format.has_key("duplicates_key") and format["duplicates_key"]:
            d = utils.removeDuplicatesFromListOfDicts(temp_data, format["duplicates_key"])
        else:
            d = temp_data
        self.linearData = d # preserves the order of the data.
        self._optimiseData()
        if not config.SILENT: print "Info: Loaded %s file: Found %s elements" % (filename, len(self.linearData))
        self._history.append("loadCSV filename=%s" % filename)
        return(True)

    def _guessDataType(self, value):
        """
        (Internal)

        Take a guess at the most reasonable datatype to store value as.
        returns the resulting data type based on a list of logical cooercions
        (explaines as I fail each cooercion).
        Used internally in _loadCSV()
        I expect this wil get larger and larger with new datatypes, so it's here as
        as separate proc.

        Datatype cooercion preference:
        float
        int
        location
        string
        """
        try: # see if the element is a float()
            return(float(value))
        except:
            try: # see if it's actually an int?
                return(int(value))
            except:
                try: # see if I can cooerce it into a location:
                    return(location(value))
                except: # this is not working, just store it as a string
                    return(str(value))

    def _optimiseData(self):
        """
        (Internal)
        Call me after modifying the data to bin and build the internal tables.
        """
        self.dataByChr = None
        if not self.linearData: # list is empty, not possible to do
            return(False)

        if self.linearData[0].has_key("loc"): # just checking the first entry.
            self.dataByChr = {}
            self.dataByChrIndexLookBack = {}
            for n, item in enumerate(self.linearData): # build the chromosome quick search maps.
                chr = item["loc"]["chr"]
                if not self.dataByChr.has_key(chr):
                    self.dataByChr[chr] = []
                    self.dataByChrIndexLookBack[chr] = []
                self.dataByChr[chr].append(item)
                self.dataByChrIndexLookBack[chr].append(n) # oh sweet, sweet dirty hack...

        # other quickmaps.
        # dataByKey?
        # end;
        return(True)

    def isChromosomeAvailable(self, chromosome):
        """
        you must check me before trying to access dataByChr[]
        """
        if self.dataByChr.has_key(chromosome):
            return(True)
        else:
            return(False)
        return(False)

    def getAllUnorderedData(self):
        """
        (Obselete?)
        return the list of chr ordered data
        """
        return(self.linearData)

    def _findDataByKeyGreedy(self, key, value): # override????? surely find?
        if key == "loc":
            # special case for locations;
            pass
        else:
            for item in self.linearData:
                if item[key] == value:
                    return(item)
        return(False)

    def _findDataByKeySlow(self, key, value): # override????? surely finditer?
        """
        finds all - returns a list
        """
        ret = []
        if key == "loc":
            # special case for locations;
            pass
        else:
            for item in self.linearData:
                if item[key] == value:
                    ret.append(item)
        return(ret)

    def getKeys(self):
        """
        return a list of all the valid keys for this geneList
        """
        return([key for key in self.linearData[0]])

    def _findAllLabelsByKey(self, key):
        """
        (Internal)
        Returns a 1D list of all the labels under Key.
        Most useful for things like geneList["Symbol"]
        geneList["entrez"]
        """
        return([x[key] for x in self.linearData])

    def _findByLabel(self, key, toFind):
        """
        (Internal)
        key is the key to search with
        toFind is some sort value to compare
        """
        for line in self.linearData:
            try:
                if toFind == line[key]:
                    return(line)
            except ValueError:
                pass
        return(None) # not found;

    def _findByCoords(self, key, coords):
        """
        (internal)
        coords should be in formal format, chrX:int(left)-int(right)
        key is unused
        """
        if not coords: return(None)

        ret = []

        for line in self.data[coords["chr"]]:
            if utils.qcollide(coords["left"], coords["right"], line["left"], line["right"]):
                ret.append(line)
        return(ret)

    def saveCSV(self, filename, **kargs):
        """
        save the geneList as a csv
        Note: This is not always available.
        As the geneList becomes more complex it loses the ability to be
        expressed simply as a csv-file. In those cases you must use
        the save() method to save a binary representation.

        **Arguments**

        filename
            filename to save, including the full path.

        **Result**

        returns None.
        saves a CSV representation of the geneList.
        """
        oh = open(filename, "w")
        writer = csv.writer(oh)
        toWrite = self.getKeys()
        writer.writerow(toWrite) # write the header row.

        for data in self.linearData:
            line = []
            for key in toWrite: # this preserves the order of the dict.
                if data.has_key(key):
                    line.append(data[key])
                else:
                    line.append("") # a blank key, fail gracefully.
            writer.writerow(line)
        oh.close()
        return(None)

    def saveFASTA(self, **kargs):
        """
        **Purpose**

        Save data as a FASTA file.

        **Arguments**

        filename
            absolute path to the file (including path)

        **Result**

        returns True if complete.
        Saves a fasta file of the sequence data in this list.

        """

        assert self.linearData[0]["seq"], "No sequence data available in this list"

        oh = open(kargs["filename"], "w")

        for item in self:
            if self.getKeys():
                name = "".join([str(item[n]) for n in self.getKeys() if not n == "seq"])
            else:
                name = "dummy"
            oh.write(">%s\r\n" % name)
            oh.write("%s\r\n" % item["seq"])
        oh.close()
        print "Info: Saved a FASTA file: %s" % kargs["filename"]
        return(True)

    def save(self, filename, **kargs):
        """
        **Purpose**

        Save the geneList as a binary representation.
        This is guaranteed to be available for all geneList representations.
        this method is used in caching the file.
        use list = load("path/to/filename")

        **Arguments**

        filename
            path to the file (including path)

        compressed
            use compression (not currently implemented)

        **Result**

        returns None
        Saves a binary representation of the geneList

        """
        compressed = False
        if kargs.has_key("compressed"): compressed = kargs["compressed"]

        oh = open(filename, "wb")
        if compressed:
            cPickle.dump(self, oh, -1)
        else:
            cPickle.dump(self, oh, -1)
        oh.close()
        print "Info: Saved Binary version of list: %s" % filename

    def sort(self, key=""):
        """
        Sort the data into a particular order based on key.

        **Arguments**

        key
            the key to use to sort by.
            must be some sort of sortable item

        **Result**

        returns True if it completes.
        sorts the list IN PLACE.
        """
        self._history.append("sort, key=%s" % key)

        assert key, "No such key"
        assert self.linearData[0].has_key(key), "Data does not have this key: %s" % key

        self.linearData = sorted(self.linearData, key=itemgetter(key))
        self._optimiseData()
        return(True)

    def reverse(self):
        """
        reverse the order of the list, in place.

        **Arguments**

        None

        **Result**

        returns True if okay or false.
        """
        self._history.append("List reversed" )
        self.linearData.reverse()
        self._optimiseData() # just in case.
        return(True)

    def getValuesInRange(self, key=None, low=None, high=None):
        """
        make a new list with all values in key that are between low and high
        """
        assert key, "must specify a key from which to extract values"
        assert low, "'low' number not valid"
        assert high, "'high' number not valid"

        newl = self.__copy__()
        newl.linearData = []
        for datum in self.linearData:
            if datum[key] > low and datum[key] < high:
                newl.linearData.append(datum)
        newl._optimiseData()
        return(newl)

    #------------------------------ Overrides --------------------------

    def __len__(self):
        """
        (Override)
        get the length of the list
        """
        return(len(self.linearData))

    def __int__(self):
        """
        (Override)
        get the length of the list
        """
        return(len(self.linearData))

    def __iter__(self):
        """
        (Override)
        make the geneList behave like a normal iterator (list)
        """
        try:
            for n in self.linearData:
                yield n
        except:
            raise StopIteration

    def __str__(self):
        """
        (Override)
        give a sensible print out.
        """
        if len(self.linearData) > config.NUM_ITEMS_TO_PRINT:
            out = []
            # welcome to perl
            for index in xrange(config.NUM_ITEMS_TO_PRINT):
                out.append("%s: %s" % (index, ", ".join(["%s: %s" % (key, self.linearData[index][key]) for key in self.linearData[index]])))
            out = "%s\r\n... truncated, showing %s/%s" % ("\r\n".join(out), config.NUM_ITEMS_TO_PRINT, len(self.linearData)+1)

            if config.PRINT_LAST_ITEM:
                out = "%s\r\n%s" % (out, "%s: %s" % (len(self.linearData), ", ".join(["%s: %s" % (key, self.linearData[-1][key]) for key in self.linearData[-1]])))

            return(out)
        elif len(self.linearData) == 0:
            return("This list is Empty")
        else: # just print first entry.
            out = []
            for index in xrange(len(self.linearData)):
                out.append("%s: %s" % (index, ", ".join(["%s: %s" % (key, self.linearData[index][key]) for key in self.linearData[index]])))
            out = "%s\r\nShowing %s/%s" % ("\r\n".join(out), len(self.linearData), len(self.linearData)+1)

            return(out)

    def __repr__(self):
        """
        (Override)
        report the underlying representation
        """
        return("glbase.genelist")

    def __getitem__(self, index):
        """
        (Override)
        confers a = geneList[0] behaviour

        This is a very slow way to access the data...
        """
        try:
            newl = self.__copy__()
            newl.linearData = []
            newl.linearData = self.linearData[index]
            newl._optimiseData()
            return(newl)
        except:
            # it is a string?
            if isinstance(index, str):
                return(self._findAllLabelsByKey(index))
            elif isinstance(index, int):
                return(self.linearData[index])
        return(False)

    def __setitem__(self, index):
        """
        (Override)
        Block key editing.
        """
        raise AssertionError, "Cannot modify list in-place"

    def __and__(self, gene_list):
        """
        (Override)
        confer and like behaviour: c = a & b
        """
        if not self.__eq__(gene_list): return(geneList()) # returns an empty list.
        newl = self.__copy__()
        newl.linearData = []
        for item1 in self.linearData:
            for item2 in gene_list.linearData:
                if item1 == item2:
                    newl.linearData.append(copy.deepcopy(item1))

        newl._optimiseData()
        self._history.append("Logical and %s & %s" % (self.name, gene_list.name))
        return(newl)

    def __or__(self, gene_list):
        """
        (Override)
        confer append like behaviour: c = a | b
        OR does not keep duplicates.
        """
        if not self.__eq__(gene_list): return(geneList())
        newl = self.__copy__()
        alist = self.linearData + gene_list.linearData
        # remove conserved duplicates;
        ulist = []
        newl = self.__copy__()
        for item in alist:
            if item in ulist:
                pass
            else:
                ulist.append(item)
                newl.linearData.append(copy.deepcopy(item))
        newl._optimiseData()
        self._history.append("Logical OR %s | %s" % (self.name, gene_list.name))
        return(newl)

    def _collectIdenticalKeys(self, gene_list):
        """
        (Internal)
        What it says, returns a list of valid keys in common between this list and gene_list
        """
        matchingKeys = []
        for key in self.linearData[0]:
            if key in gene_list.linearData[0]:
                matchingKeys.append(key)
        return(matchingKeys)

    def __add__(self, gene_list):
        """
        (Override)
        confer append like behaviour: c = a + b
        keeps duplicates (just concatenate's lists)
        """
        mkeys = self._collectIdenticalKeys(gene_list)
        if not mkeys: # unable to match.
            print "Warning: No matching keys, the resulting list would be meaningless"
            return(False)
        newl = self.__copy__()
        newl.linearData = []
        newl.linearData = self.linearData + copy.deepcopy(gene_list.linearData)
        newl._optimiseData()
        self._history.append("List addition %s + %s" % (self.name, gene_list.name))
        return(newl)

    def __sub__(self, gene_list):
        """
        (Override)
        confer c = a - b ability.
        Actually xor?
        """
        mkeys = self._collectIdenticalKeys(gene_list)
        if not mkeys: # unable to match.
            print "Warning: No matching keys, unable to perform subtraction"
            return(False)

        newl = self.__copy__()
        newl.linearData = []

        dontAdd = False
        for item in self.linearData: # do a map here...
            for item2 in gene_list.linearData:
                for k in mkeys:
                    if item[k] == item2[k]:
                        dontAdd = True
                    else:
                        dontAdd = False # all mkeys must match
            if not dontAdd:
                newl.linearData.append(copy.deepcopy(item))
            dontAdd = False
        newl._optimiseData()
        self._history.append("List subtraction %s - %s" % (self.name, gene_list.name))
        return(newl)

    def __eq__(self, gene_list):
        """
        (Internal)
        Are the lists equivalent?
        lists now, must only have one identical key.
        """
        # check the hash's first to see if they are identical.
        if self.__hash__() == gene_list.__hash__():
            return(True)

        for key in self.linearData[0]:
            if key in gene_list.linearData[0]:
                return(True) # just one key in common required.
        if not config.SILENT: print "Warning: gene lists are not equivalent, (you are comparing non-identical format lists)"
        return(False)

    def __hash__(self):
        """
        (Override)

        compute a sensible hash value
        """
        try:
            return(hash(self.name + str(self[0]) + str(len(self)))) # hash data for comparison.
        except:
            try:
                return(hash(self.name + str(self[0]))) # len() probably not available.
            except: # I bet the list is empty.
                return(hash(self.name))


    def __ne__(self, gene_list):
        """
        (Internal)
        Are the lists equivalent?
        ie do they have the same keys?
        """
        return(not self.__eq__(gene_list))

    def __copy__(self):
        """
        (Override)
        Confer copy to mean a deep copy as opposed to a shallow copy.
        """
        return(copy.deepcopy(self))

    def append(self, item=None, **kargs):
        """
        Append an item onto the geneList.
        The appended item must have all valid keys to work.
        To test if the item can be appended to this list call::

            if list == item: # test they are compatible.
                list.append(item) # append to list.

        NOTE: Be careful in it's usage, it optimises the data on each
        append, potentially slow with large numbers of appends()

        In these cases it is better to add the lists together::

            new_list = old_list + other_list # A simple append, no duplicate removal
            new_list = old_list & other_list # AND the lists together (only keep items in both lists)
            new_list = old_list | other_list # OR the lists, keep all items that appear in either list.

        **Arguments**

        item
            The item to append. It must be the same format as the list
            to append()

        **Result**

        returns None

        Updates the list in-place.
        """
        assert item, "Cannot add an empty item onto the list"
        # Don't I want to test for geneList like structures?
        # Test for equivalence of keys?
        self.linearData.append(item) # this just appends to the list?
        self._optimiseData() # oof!!!

    def getColumns(self, return_keys=None):
        """
        return a new geneList only containing the coluns specified in return _keys (a list)
        """
        if not isinstance(return_keys, list):
            return_keys = []
            return_keys.append(return_keys)

        newl = self.__copy__()
        newl.linearData = []
        for item in self.linearData:
            newd = {}
            for key in return_keys:
                newd[key] = item[key]
            newl.linearData.append(newd)
        newl._optimiseData()
        self._history.append("Column sliced, for the columns: %s" % return_keys)
        return(newl)

    def getCategories(self):
        """
        Get all of the categories in the list;
        """
        return([key for key in self.linearData[0]])

    def map(self, genelist=None, peaklist=None, microarray=None, genome=None, key=None, **kargs):
        """
        **Purpose**

        map() is an important method, it merges merges two gene_list-like
        objects and outputs a new genelist.

        The new gene_list will inherit from 'the right', for
        example if you have a microarray you should perform the map in this
        order::

            result = gene_list.map(microarray=microarray, "refseq")

        'result' will now be a microarray object with all the appropriate
        methods.

        If however, you did this by mistake::

            result = microarray.map(genelist=gene_list, "refseq")

        It will still work fine, but now, trying::

            result.drawHeatmap(...)

        will fail, because the result is a vanilla gene_list rather than
        a microarray as you might intend.

        Also note, the algortihm is 'greedy' and will only take the first
        matching entry it finds.

        **Arguments**

        genelist or peaklist or microarray
            some sort of inherited gene_list like object,
            examples inglude genelist, microarray, genome, taglist
            note: seqfile is not a genelist object and cannot
            be used in this manner.

        key
            a key in common between the two lists you can use to map
            them against each other.

        **Result**

        returns a new gene_list-like object, inheriting methods from the
        left hand side of the equation.

        """
        assert genome or genelist or peaklist or microarray, "No valid genelist specified"
        if genelist: gene_list = genelist
        if peaklist: gene_list = peaklist
        if microarray: gene_list = microarray
        if genome: gene_list = genome
        assert key, "Must specify a 'key' to map the two lists"
        assert gene_list.linearData[0].has_key(key), "the genelist must have the key"
        map_key = key

        # speed up with a triangular search?
        newl = gene_list.__copy__()
        newl.linearData = []
        for item in self.linearData:
            for other in gene_list:
                #print item, other
                if item[map_key] == other[map_key]:
                    new_entry = copy.deepcopy(item)
                    for key in other:
                        if key != map_key:
                            new_entry[key] = other[key]
                    newl.linearData.append(new_entry)
                    break
        if not config.SILENT: print "Info: Mapped two lists by %s, found: %s" % (map_key, len(newl))
        if len(newl):
            newl._optimiseData()
            return(newl)
        else:
            if not config.SILENT: print "Warning: No keys matched using: %s" % map_key
            return(False)

    def _findNearbyGenes(self, coords, distance=10000):
        """
        expects:
        coords = chr1:10000-10002
        distance = distance from the coords;
        # relies on you list having a valid tss_loc
        """
        assert coords, "Cannot annotate: %s" % coords
        assert self.linearData[0].has_key("tss_loc"), "Error: no available tss_loc key"

        if not self.isChromosomeAvailable(str(coords["chr"])): return(False) # has this chromsome.

        peakCentre = (coords["left"] + coords["right"]) / 2

        ret = []
        for line in self.dataByChr[str(coords["chr"])]:
            line_loc = line["tss_loc"]
            tss_start = line_loc["left"]

            if (peakCentre >= (tss_start-distance)) and (peakCentre <= (tss_start+distance)):
                # add a new _dist_to_tss tag;
                line["dist_to_tss"] = peakCentre-tss_start
                #print line, peakCentre, tss_start
                ret.append(line)
        return(ret)

    def annotate(self, gene_list=None, key_to_match="loc", distance=10000, resolution=2000, **kargs):
        """
        Annotate this gene_list by referring to some other gene_list-like
        object, often a genome. It will map loc within a certain distance

        **Arguments**

        gene_list
            A genelist-like object

        key_to_match
            must be some sort of location tag::

            e.g.
            * loc "chr1:1110100-1110200"
            * location "chr1:1110100-1110200"
            * tss_loc "chr1:1110100-1110200"
            * etc.

        distance
            The distance (in base pairs) to look for a match.

        resolution (Optional)
            annotate() draws an image of the moving window and saves it
            in the "gene_list"'s path.

            This image will have the filename "{list_name}_tss_dist_distribution.png"
            the resolution argument specifies the size of the moving window
            used in the calculation of the graph.

        image_filename (Optional)


        **Result**

        returns a new gene_list-like object inherited from the right-side.
        i.e.::

            result = gene_list.map(microarry, "tss_loc")

        result will be a microarray.

        however::

            result = microarry.map(chip_list, "loc")

        result will be a chip_list

        Saves the frequency of hits relative to the location tag as a PNG image
        in the current path with the filename: "{list_name}_tss_dist_distribution.png"
        """
        # -------------- Error Checking
        if not config.SILENT: print "Info: Annotate %s..." % gene_list.name

        assert gene_list, "you must specify a gene_list"
        assert gene_list.linearData[0].has_key(key_to_match), "Cannot annotate, no tss_location, do you have the genome and gene list the wrong way around?"
        assert self.linearData[0].has_key("tss_loc"), "The annotation list does not have a valid transcription start-site key"

        # --------------- Parse the kargs
        image_filename = None
        if kargs.has_key("image_filename"): png_filename = kargs["png_filename"]

        # --------------- Proper
        newl = gene_list.__copy__() # make a new copy, I need to copy to preserve the attributes.
        # or if it is a derived class.
        newl.linearData = [] # blank the data though

        headerLabels = set(["chr", "loc", "name"])
        dist_hist = []

        total_hits = 0

        for index, entry in enumerate(gene_list):
            hits = self._findNearbyGenes(entry[key_to_match], distance)
            if hits:
                total_hits += len(hits)
                if not config.SILENT and config.VERBOSE: print "Found: %s hits for entry: %s/%s" % (len(hits), index, len(gene_list))
                for annotation in hits:
                    # I want to merge the two lists in a new dict;
                    new_entry = {}
                    for key in annotation:
                        new_entry[key] = annotation[key]
                    for key in entry:
                        new_entry[key] = entry[key]
                    newl.linearData.append(new_entry)
                    if image_filename:
                        dist_hist.append(int(annotation["dist_to_tss"]))

        if image_filename:
            window = resolution
            # first I need to bin the data (by 1) - then do the moving average.
            linData = [0 for x in xrange(-distance, distance)] # set-up a blank array
            for item in dist_hist:
                linData[item-distance] += 1
            x, m_average = utils.movingAverage(linData, window)

            self.draw._plot(os.path.join(gene_list.path, png_filename),
                m_average, x=x,
                title="%s - loc of Chipseq around tss, moving window = %s" % (gene_list.name, window),
                xlabel="Distance to transcription start site (base pairs)",
                ylabel="Frequency (Raw)")

        if not config.SILENT: print "Info: Annotated %s, found: %s" % (newl.name, len(newl))
        newl._optimiseData()
        self._history.append("Annotate list, %s using: %s, using key: %s" % (self.name, gene_list.name, key_to_match))
        return(newl)

    def collide(self, loc_key="loc", delta=200, **kargs):
        """
        **Purpose**

        collide two lists using some sort of location tag or key.
        takes two lists, the parent list and one of microarray, peaklist
        or other geneList-like object with a "location"-like key.

        Be careful in your usage of collide() and overlap(), they both merge
        coordinates, but collide() uses the centre of the peak, whereas
        overlap() uses the centre's of the peaks.

        **Arguments**

        genelist or peaklist or microarray
            must be some sort of geneList-like object

        loc_key (Optional, Default: "loc")
            the key to use as a location tag.

        delta (Optional, Default: 200)
            the collision fuzzy window +- around the middle of the loc_key tag

        logic (Optional, Default: "and")
            can be one of "and", "not"
            "and" = keep only collisions in both lists.
            "not" = perform a 'not collision', keeping elements that do not collide.
            "notleft" = keep only the left side (i.e. the parent list)
            "notright" = keep only the right list (i.e. peaklist|genelist|microarray argument"

        image_filename (Optional)
            save a dot-plot map of the colliding lists.

        merge (Optional, Default: False)
            merge the two lists keys together. This is useful for
            heterogenouse lists. For homogenous lists, don't set this
            or set it to False and collide() will keep the resulting list
            homogenous too.
            This will only merge 'extra' keys, it does not overwrite the
            original lists keys where two values exist.

        **Result**

        returns a new genelist derived from the intersect of the keys in
        the original list and the new gene list
        the resulting "loc" key will be a single base pair mid-way
        between the centres of the two colliding coordinates.
        """
        if kargs.has_key("peaklist"): gene_list = kargs["peaklist"]
        elif kargs.has_key("genelist"): gene_list = kargs["genelist"]
        elif kargs.has_key("microarray"): gene_list = kargs["microarray"]
        else:
            print "Error: No valid gene_list"
            quit()

        assert self[0].has_key(loc_key), "Error: no location key"
        assert gene_list[0].has_key(loc_key), "Error: no location key"

        merge = False
        if kargs.has_key("merge"):
            merge = True

        mode = "and"
        if kargs.has_key("logic") and kargs["logic"] != "and":
            mode = kargs["logic"]
            foundA = [False for x in xrange(len(self))]
            foundB = [False for x in xrange(len(gene_list))]

        newl = gene_list.__copy__()
        newl.linearData = []
        if kargs.has_key("image_filename"):
            raise NotImplementedError, "Plotting the dotplot of collide is not implemented yet"
            result_array = [0 for x in xrange(max(len(gene_list), len(self)))]
            x_data = []
            y_data = []

        for indexA, item in enumerate(self):
            locA = item[loc_key]
            if gene_list.isChromosomeAvailable(locA["chr"]):
                Acentre = (locA["right"] + locA["left"]) / 2

                for local_index, other in enumerate(gene_list.dataByChr[locA["chr"]]):
                    indexB = gene_list.dataByChrIndexLookBack[locA["chr"]][local_index] # sweet sweet hack, see _optimiseData for where this comes from
                    locB = other[loc_key]
                    Bcentre = (locB["right"] + locB["left"]) / 2

                    if utils.collide(Acentre - delta, Acentre + delta, Bcentre - delta, Bcentre + delta):
                        if mode != "and":
                            foundB[indexB] = True # only matters if logic="not"
                            foundA[indexA] = True
                        if merge:
                            a = copy.deepcopy(item)
                            b = copy.deepcopy(other)
                            for k in b:
                                if not a.has_key(k):
                                    a[k] = b[k]
                        else:
                            a = copy.deepcopy(item)

                        a["dist"] = Acentre - Bcentre # add a new key.
                        a["loc"] = location(loc=locA)
                        a["loc"].pointify()

                        for tag_key in ["tags", "tag", "tag_height"]: # if I find these keys I will add them together.
                            if item.has_key(tag_key) and other.has_key(tag_key):
                                try:
                                    a[tag_key] = int(item[tag_key]) + int(other[tag_key])
                                except:
                                    a[tag_key] = item[tag_key]

                        newl.linearData.append(a)
                        if kargs.has_key("png_filename"):
                            result_array[indexA] += 1
                            result_array[indexB] += 1
                            x_data.append(indexA)
                            y_data.append(indexB)

        # if logic = "not", I want the opposite of all the items found
        if mode in ["not", "notleft", "notright"]:
            newl.linearData = []
            if mode == "not" or mode == "notleft":
                for index, item in enumerate(self):
                    if not foundA[index]:
                        a = copy.deepcopy(item)
                        loc = a["loc"]
                        a["loc"] = location(loc=loc) # pointify result
                        a["loc"].pointify()
                        newl.linearData.append(a)

            if mode == "not" or mode == "notright":
                for index, item in enumerate(gene_list):
                    if not foundB[index]:
                        a = copy.deepcopy(item)
                        loc = a["loc"]
                        a["loc"] = location(loc=loc) # pointify result
                        a["loc"].pointify()
                        newl.linearData.append(a)
        """
        if karg.has_key("image_filename"):
            x, data = utils.movingAverage(result_array, 100)
            plot.cla()
            plot.plot(data)
            plot.savefig("meh.png")

            plot.cla()
            plot.scatter(x_data, y_data)
            plot.savefig("sca.png")
        """
        if not config.SILENT:
            if mode == "and":
                 print "Info: Collided two lists using %s, found: %s overlaps" % (loc_key, len(newl))
            elif mode == "not":
                print "Info: \"not\" Collided two lists using %s, found: %s non-overlapping sites" % (loc_key, len(newl))
        newl._optimiseData()
        self._history.append("List Collision %s, using: %s, with key: %s" % (self.name, gene_list.name, loc_key))
        return(newl)

    def overlap(self, loc_key="loc", delta=200, **kargs):
        """
        **Purpose**

        overlap two lists using some sort of location tag or key.
        takes two lists, the parent list and one of microarray, peaklist
        or other geneList-like object with a "location"-like key.

        Be careful in your usage of collide() and overlap(), they both merge
        coordinates, but collide() uses the centre of the peak, whereas
        overlap() uses the centre's of the peaks.

        **Arguments**

        genelist or peaklist or microarray
            must be some sort of geneList-like object

        loc_key (Optional, Default: "loc")
            the key to use as a location tag.

        delta (Optional, Default: 200)
            the collision fuzzy window +- around the left and right
            points of the loc_key tag

        logic (Optional, Default: "and")
            can be one of "and", "not"
            "and" = keep only collisions in both lists.
            "not" = perform a 'not collision', keeping elements that do not collide.
            "notleft" = keep only the left side (i.e. the parent list)
            "notright" = keep only the right list (i.e. peaklist|genelist|microarray argument"

        image_filename (Optional)
            save a dot-plot map of the colliding lists.

        merge (Optional, Default: False)
            merge the two lists keys together. This is useful for
            heterogenouse lists. For homogenous lists, don't set this
            or set it to False and collide() will keep the resulting list
            homogenous too.

        **Result**

        returns a new genelist derived from the intersect of the keys in
        the original list and the new gene list
        the resulting "loc" key will be the left and right most extreme points
        of the two overlapping coordinates.
        """
        if kargs.has_key("peaklist"): gene_list = kargs["peaklist"]
        elif kargs.has_key("genelist"): gene_list = kargs["genelist"]
        elif kargs.has_key("microarray"): gene_list = kargs["microarray"]
        else:
            print "Error: No valid gene_list"
            quit()

        assert self[0].has_key(loc_key), "Error: no location key"
        assert gene_list[0].has_key(loc_key), "Error: no location key"

        merge = False
        if kargs.has_key("merge"):
            merge = True

        mode = "and"
        if kargs.has_key("logic") and kargs["logic"] != "and":
            mode = kargs["logic"]
            foundA = [False for x in xrange(len(self))]
            foundB = [False for x in xrange(len(gene_list))]

        newl = gene_list.__copy__()
        newl.linearData = []
        if kargs.has_key("image_filename"):
            raise NotImplementeError, "Plotting the dotplot of overlap is not implemented yet"
            result_array = [0 for x in xrange(max(len(gene_list), len(self)))]
            x_data = []
            y_data = []

        for indexA, item in enumerate(self):
            locA = item[loc_key]
            if gene_list.isChromosomeAvailable(locA["chr"]):
                for local_index, other in enumerate(gene_list.dataByChr[locA["chr"]]):
                    indexB = gene_list.dataByChrIndexLookBack[locA["chr"]][local_index] # sweet sweet hack, see _optimiseData for where this comes from
                    locB = other[loc_key]

                    if utils.collide(locA["left"]-delta, locA["right"]+delta, locB["left"]-delta, locB["right"]+delta):
                        if mode != "and":
                            foundB[indexB] = True # only matters if logic="not"
                            foundA[indexA] = True
                        if merge:
                            a = copy.deepcopy(item)
                            b = copy.deepcopy(other)
                            for k in b:
                                if not a.has_key(k):
                                    a[k] = b[k]
                        else:
                            a = copy.deepcopy(item)

                        a["loc"] = location(chr=locA["chr"], left=min(locA["left"], locB["left"]), right=max(locA["right"], locB["right"]))

                        for tag_key in ["tags", "tag", "tag_height"]: # if I find these keys I will add them together.
                            if item.has_key(tag_key) and other.has_key(tag_key):
                                try:
                                    a[tag_key] = int(item[tag_key]) + int(other[tag_key])
                                except:
                                    a[tag_key] = item[tag_key]

                        newl.linearData.append(a)
                        if kargs.has_key("png_filename"):
                            result_array[indexA] += 1
                            result_array[indexB] += 1
                            x_data.append(indexA)
                            y_data.append(indexB)

        # if logic = "not", I want the opposite of all the items found
        if mode in ["not", "notleft", "notright"]:
            newl.linearData = []
            if mode == "not" or mode == "notleft":
                for index, item in enumerate(self):
                    if not foundA[index]:
                        a = copy.deepcopy(item)
                        a["loc"] = location(loc = a["loc"]) # pointify result
                        a["loc"].pointify()
                        newl.linearData.append(a)

            if mode == "not" or mode == "notright":
                for index, item in enumerate(gene_list):
                    if not foundB[index]:
                        a = copy.deepcopy(item)
                        a["loc"] = location(loc=a["loc"]) # pointify result
                        a["loc"].pointify()
                        newl.linearData.append(a)
        """
        if karg.has_key("image_filename"):
            x, data = utils.movingAverage(result_array, 100)
            plot.cla()
            plot.plot(data)
            plot.savefig("meh.png")

            plot.cla()
            plot.scatter(x_data, y_data)
            plot.savefig("sca.png")
        """
        if not config.SILENT:
            if mode == "and":
                 print "Info: Collided two lists using %s, found: %s overlaps" % (loc_key, len(newl))
            elif mode == "not":
                print "Info: \"not\" Collided two lists using %s, found: %s non-overlapping sites" % (loc_key, len(newl))
        newl._optimiseData()
        self._history.append("List Collision %s, using: %s, with key: %s" % (self.name, gene_list.name, loc_key))
        return(newl)

    def history(self):
        """
        get the origins and history and where this list has been, and all of the operators
        performed on it.
        """
        print self._history

    def removeDuplicates(self, key=None):
        """
        remove the duplicates in the list and returns a new list;
        only keeps the first example (greedy)
        """
        assert key, "No key specified"
        newl = self.__copy__()
        newl.linearData = []
        ulist = []
        for item in self:
            if not item[key] in ulist:
                ulist.append(item[key])
                newl.linearData.append(copy.deepcopy(item))
        newl._optimiseData()
        return(newl)
