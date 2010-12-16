"""
Genome
======

TODO:
-----
* sequence should be bindable from the __init__()

"""

import sys, os, csv
import config, utils

from genelist import genelist
from flags import *
from errors import AssertionError
from location import location

class genome(genelist):
    def __init__(self, name="None", filename=None, **kargs):
        """
        **Purpose**

        genome container.

        **Arguments**

        filename
            absolute filename (including path) to the actual file.
            can include path short cuts (e.g. "./", "../" etc)

        format (Optional, default = "sniffer" (ie. guess))
            format specifer, see flags.py and helpers.py and the
            documentation on how to write a valid format specifier

        sequence_path (Optional, default=None)
            provide a path to the fasta sequence files.
            If you provide the fasta chromosome files for the
            correct genome assembly then you can extract genomic
            sequence using getSequence() and getSequences()

        **Result**

        fills the genelist with the data from the file as specified by
        the format specifier.
        If no filename is specified it returns an empty genome.
        """
        valig_args = ["filename", "format", "sequence_path"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.__init__, k)

        genelist.__init__(self) # inherit

        if filename:
            self.name = name
            self.filename = filename
            if "format" in kargs:
                format = kargs["format"]
            else:
                format = default # sniffer

            # sniff the file to see if it's a tsv or a csv.
            if "sniffer" in format:
                oh = open(filename, "rU")
                reader = csv.reader(oh)
                for item in reader:
                    if len(item) == 1:
                        format["dialect"] = csv.excel_tab
                    break
                oh.close()

            self.loadCSV(filename=filename, format=format) # use the standard loader.

            if not "tss_loc" in self.linearData[0]:
                # the list does not have a tss_loc key,
                # I need to build it myself from the loc and strand tags
                assert "loc" in self.linearData[0], "I can't find a location in the genome data! you must specify a 'tss_loc' or 'loc' key"

                for item in self.linearData:
                    cloc = item["loc"]
                    if item["strand"] in positive_strand_labels:
                        item["tss_loc"] = location(chr=cloc["chr"], left=cloc["left"], right=cloc["left"])
                    elif item["strand"] in negative_strand_labels:
                        item["tss_loc"] = location(chr=cloc["chr"], left=cloc["right"], right=cloc["right"])

            self._optimiseData()
        if "sequence_path" in kargs:
            self.bindSequence(kargs["sequence_path"])
        else:
            self.bHasBoundSequence = False

        self._optimiseData()

    def __str__(self):
        """
        (Internal)
        Override the genelist version to append a flag
        showing presence (or not) of bound genomic sequence.
        """
        ret = genelist.__str__(self)
        if self.bHasBoundSequence:
            return("%s\n(DNA sequence is available)" % ret)
        else:
            return("%s\n(DNA sequence is not available)" % ret)

    def __repr__(self):
        return("glbase.genome")

    def findGenes(self, **args):
        """
        returns an appended normal vanilla list
        candidate for internalisation?
        valid args:
        loc = get gene by coords, they must be exact, returns the formal coords and the refseqID
        refseq = get the gene info by RefSeqID, returns the gene coords
        key = find by the value key...

        although you could pass both at the same time, this will only respect the first arg.
        make it respectable?
        """
        #print "a",args
        ret = []

        func_dict = {
            "refseq": self._findByLabel, # these procedures will always accept two vars, they may not use them though.
            "loc": self._findByCoords,
            "tss_loc": self._findByCoords,
            "entrez": self._findByLabel
            }
        for key, value in args.iteritems():
            if key == "key": # not a documented way to do it, but if you want to just search by key. pass a tuple of the form (actual_key, data)
                if value[0] in self.linearData[0]:
                    return(func_dict[value[0]](value[0], value[1]))
                else:
                    return(None)
            else: # the normal documented way
                if key in self.linearData[0]:
                    func_dict[key](key, value)
            return(None)

    def getAnnotationsByKey(self, geneList, key, KeepAll=False):
        """
        get a list of annotations from a geneList based on key that is common between the two lists.
        keepAll == keep even the unmapped ones.
        """
        assert key in self.linearData[0], "Error: Genome does not have search key: %s" % key
        assert key in geneList.linearData[0], "Error: Genome does not have search key: %s" % key

        omit = 0

        mockEntry = {}
        for key in self.linearData[0]: # make a mock entry for the KeepAll Situation.
            mockEntry[key] = ""

        newl = []
        for item in geneList.linearData:
            gene = self.findGenes(key=(key, item[key]))
            if gene:
                newl.append(gene)
            else:
                omit += 1
                if KeepUnmapped: # send back a fake return
                    newEmptyEntry = mockEntry
                    newEmptyEntry[key] = item[key]
                    newl.append(newEmptyEntry)
        config.log.error("Could not find: %s/%s" % (omit, len(geneList)))
        return(newl)

    def getDataByFeature(self, geneList, feature="refseq", returnFeatures=["name", "refseq", "entrez", "tss_loc"], KeepUnmapped=False):
        """
        returns the coords of genes based on geneList and returns a list of features.
        if geneList=None, it will return all genes, (you can use this to extract sets)
        This is different from getAnnotationsByKey in that it only returns a list in the form returnFeatures.
        This makes it convenient to dump into a csv file.
        """
        if not geneList:
            ret = []
            for item in self.linearData:
                ret.append([item[feature], item["name"]])
            return(ret)

        ret = self.getAnnotationsByKey(geneList, feature, KeepUnmapped) # get all the items

        new_ret = []
        if ret:
            for item in ret:
                row = []
                for key in returnFeatures:
                    row.append(item[key])
                new_ret.append(row)
        return(new_ret)

    def getFeatures(self, loc=None, **kargs):
        """
        **Purpose**
            get all of the genomic features (probably genes) within a
            certain location span.

        **Arguments**
            location
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
                    item["type"] = "gene" # set the type flag for gDraw
                    ret.append(item)
        return(ret)

    def bindSequence(self, path=None):
        """
        **Purpose**

        Bind genome fasta files so that this genome object will recognise
        the sequence. This step is required if you want to use fastalists
        and genome.getSequence()

        **Arguments**

        path
            path specifying the locations of the FASTA files that make
            up the sequence data. They usually come in the form "chr1.fa"
            for human and mouse genomes.

        **Result**

        returns True if complete.
        genome.getSequence(loc="chrN:left-right") will now correctly return
        the sequence specified by the location.
        """
        assert os.path.realpath(path), "'path' argument is required"
        assert os.path.exists(os.path.realpath(path)), "Path does not exist"

        self.fasta_dir = os.path.realpath(path)
        fastas = os.listdir(self.fasta_dir)
        assert len(fastas) > 0, "no fasta files found"
        self.seq = {}
        self.seq_data = {}
        for f in fastas:
            if f[0:3] == "chr" and f[-3:] == ".fa":
                chr_key = f.split(".")[0].replace("chr", "")
                self.seq[chr_key] = open(os.path.join(path, f), "rU")
                self.seq_data[chr_key] = {"filename": f, "offset": 0, "linelength": 0}

                # open the file and sniff the fasta data.
                first_line = self.seq[chr_key].readline() # should be a fasta handle >dddd d
                assert first_line[0] == ">", "Not a valid FASTA file, header is missing"
                self.seq_data[chr_key]["offset"] = self.seq[chr_key].tell() # get the location after the fasta header
                self.seq_data[chr_key]["fasta_header"] = first_line

                second_line = self.seq[chr_key].readline() # should be sequence.
                assert len(second_line[0]) > 0, "Not a valid FASTA file, sequence is missing?"
                self.seq_data[chr_key]["linelength"] = len(second_line) -1 # each line.
        config.log.info("Bound Genome sequence: %s" % path)
        self.bHasBoundSequence = True
        return(True)

    def getSequence(self, **kargs):
        """
        **Purpose**

        get the sequence under coords...

        **Arguments**

        coords or loc (one is Required)
            genomic coordinates of the form "chr1:100100-100200"

        strand
            Use, +, f, top, forward, 0, for the top strand
            Use, -, r, bottom, reverse, 1, for the reverse complement strand
            If the strand is not specified then the + strand will be returned.

        mask (Optional, default=False)
            'repeat mask' the returned sequence (i.e. convert lower-case
            acgt to NNNN)

        **Result**

        returns a string containing the sequence at 'coords'
        """
        valid_args = ["coords", "loc", "strand", "mask"]
        for key in kargs:
            assert key in valid_args, "getSequence() - Argument '%s' is not recognised" % key

        assert "loc" in kargs or "coords" in kargs, "No valid coords or loc specified"
        assert self.bHasBoundSequence, "No Available genome FASTA files"

        if "loc" in kargs: loc = kargs["loc"]
        elif "coords" in kargs: loc = kargs["coords"]
        try:
            loc = location(loc=loc)
        except:
            pass
        assert isinstance(loc, location), "'loc' must be a proper genome location"

        left = loc["left"]
        right = loc["right"]
        chrom = loc["chr"]

        seekloc = (left + (left / self.seq_data[chrom]["linelength"]))-1 # the division by 50 is due to the presence of newlines every 50 characters.
        self.seq[chrom].seek(seekloc+self.seq_data[chrom]["offset"]) # move to the start location.

        delta = (right - left)+1

        # count the number of line endings.
        # get a niave reading.
        bonus = 0
        ret = ""
        while len(ret) < delta:
            self.seq[chrom].seek(seekloc+self.seq_data[chrom]["offset"])
            ret = self.seq[chrom].read(delta + (delta /self.seq_data[chrom]["linelength"]) + bonus).replace("\n", "").replace("\r", "")
            bonus += 1
            if bonus > delta: # breaks in case you send a loc that is beyond the end of the file.
                break

        if "strand" in kargs:
            if kargs["strand"] in negative_strand_labels:
                ret = utils.rc(ret)

        if "mask" in kargs and kargs["mask"]:
            ret = utils.repeat_mask(ret)

        return(ret)

    def getSequences(self, **kargs):
        """
        **Purpose**

        Get all of the sequences from a gene_list-like object (e.g. a peaklist,
        microarray, etc) with some sort of valid location key (e.g. "chr1:10000-20000")

        I've checked this extensively - if you provide the correct location
        it will send back the correct sequence. Any further errors are
        generally from wrong locations. So Heh.

        **Arguments**

        genelist (Required)
            some sort of genelist-like object

        loc_key (Required)
            the name of the location key (e.g. "loc", "tss_loc", "coords", etc..)

        Optional Arguments

        strand_key
            If you want the list to respect the strand you must tell it the name of
            a 'strand' key.

        deltaleft=n
            expand the coordinates left by n base pairs
            Will respect the orientation of the strand.
            (You must specify a 'strand_key')

        deltaright=n
            expand the coordinates rightwards by n base pairs.
            Will respect the orientation
            of the strand.
            (You must specify a 'strand_key')

        delta=n
            expand the coordinates by n, added onto the left and right

        pointify (True|False, default=False)
            turn the location into a single base pair based on the centre
            of the coordinates (best used in combination with delta to
            expand reads symmetrically, pointify will be performed before
            the expansion of the coordinates)

        mask (default=False)
            use the upper and lower case of the fasta files to 'mask'
            the sequence. This will turn acgt to NNNN.

            This is not a proper repat masker and relies on your genome
            being repeat masked. For human and mouse this is usually true,
            but for other genomes has not been tested.

        **Result**

        returns a copy of the original genelist-like object
        with the new key "seq" containing the sequence.
        Will add a new key "seq_loc" that contains the new location that
        the seq spans across.
        """
        valid_args = ["genelist", "loc_key", "strand_key", "deltaleft", "deltaright", "delta",
            "pointify", "mask"]
        for key in kargs:
            assert key in valid_args, "getSequences - Argument '%s' is not recognised" % key

        assert self.bHasBoundSequence, "No Available genome FASTA files"
        assert "genelist" in kargs, "Required argument: 'list' is missing or malformed"
        assert "loc_key" in kargs, "Required argument: 'loc_key' is missing or malformed"
        #assert "deltaleft" in kargs and "strand_key" in kargs, "You must specify a strand_key if you want to assymetrically expand the coordinates"
        #assert "deltaright" in kargs and "strand_key" in kargs, "You must specify a strand_key if you want to assymetrically expand the coordinates"

        newl = kargs["genelist"].__copy__()
        loc_key = kargs["loc_key"]

        strand_key = False
        if "strand_key" in kargs:
            strand_key = kargs["strand_key"]

        mask = False
        if "mask" in kargs and kargs["mask"]:
            mask = True

        for item in newl.linearData: # this will break delayedlists...
            newloc = item[loc_key]

            if "pointify" in kargs and kargs["pointify"]:
                newloc = newloc.pointify()

            if "delta" in kargs:
                newloc = newloc.expand(kargs["delta"])

            if "deltaleft" in kargs and kargs["deltaleft"]:
                if kargs["strand_key"] in positive_strand_labels:
                    newloc = newloc.expandLeft(kargs["deltaleft"])
                elif kargs["strand_key"] in negative_strand_labels:
                    newloc = newloc.expandRight(kargs["deltaleft"])

            if "deltaright" in kargs and kargs["deltaright"]:
                if kargs["strand_key"] in positive_strand_labels:
                    newloc = newloc.expandLeft(kargs["deltaright"])
                elif kargs["strand_key"] in negative_strand_labels:
                    newloc = newloc.expandRight(kargs["deltaright"])

            if strand_key:
                seq = self.getSequence(loc=newloc, strand=item[strand_key], mask=mask)
                item["strand"] = item[strand_key]
            else:
                seq = self.getSequence(loc=newloc, mask=mask) # defaults to + strand.
                item["strand"] = "+" # overwrite the strand to reflect the seq

            item["seq"] = seq
            item["seq_loc"] = newloc

        newl._optimiseData()
        config.log.info("Got sequences for '%s'" % self.name)
        newl._history.append("Added DNA sequence")
        return(newl)
