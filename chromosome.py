"""
chromosome, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

"""

import utils, error

chrNames = ["X", "Y", "M", "x", "y", "m"] # valid chromosome names;

class chromosome:
    def __init__(self):
        """
        chromosome constructor
        """
        self.number = 0 # chromosme number/name
        self.feature = []


    def setNumber(self, _n):
        """
        set the chromosome number (or name, if X, Y, M)
        """
        try:
            self.number = int(_n)
        except ValueError:
            if _n in chrNames:
                self.number = n
            else:
                error.error("Error: 1; Bad chromosome name", True)

    def bindFASTASeqFile(self, fastaFile):
        """
        bind a raw sequence file to the chromosome
        """
        self.fastaFile = fastaFile
        self.fastaHandle = open(fastaFile, "wb")
        # parse file to find starts and end values for grabbing sequence.

    def bindFeature(self):
        """
        bind a feature set to the chromosome
        """
        pass

    def getSequenceA(self, left, right):
        pass

    def getSequenceB(self, chr_loc):
        """
        getSequence using Chr:left-right nomenclature
        """
        loc = utils.getLocation()
        if loc["chr"] == self.number:
            pass
        else:
            return("")

    def getFeatures(self, left, right):
        """
        Extract the features and send them back in a displayable format
        from left to right of the genome.
        """
        pass

