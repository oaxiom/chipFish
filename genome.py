"""
genome, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Would this be better inheriting from feature?

"""

import csv, os, sys

import utils, opt

from feature import *

spcode = {
    "Mus Musculus": "mm",
    "mouse": "mm",
    "Mouse": "mm",
    "Homo Sapiens": "hs",
    "human": "hs",
    "Human": "hs",
    "Uberalles": "hs", #  hehe :)
    "Arabidopsis": "at",
    "Arabidopsis thalina": "at",
    "Evil Weed": "at" # oh, I'm so funny... .. um....
    }

class genome:
    def __init__(self, species, assembly):
        """
        Initialise the genome,
        species is a simple species string that corresponds to one of the species codes.
        """
        self.assembly = assembly
        self.chrl = []
        self.species = species
        self.features = []
        if spcode.has_key(species):
            self.code = spcode[species.lower()]
        else:
            self.code = "uu"
        self.path = os.path.join(opt.generic.app_path, "data/", self.code, assembly)

        # load the genome meta data, no. chr, genes, etc...
        # this will be a db query later, but for now we just use csv's
        oh = open(os.path.join(self.path, "meta_data.csv"))
        reader = csv.reader(oh)
        for line in reader:
            if line[0] == "n_chromosomes":
                self.nChromosomes = int(line[1])
            elif line[0] == "hasX":
                self.hasX = True
            elif line[0] == "hasY":
                self.hasY = True

        oh.close()

        self.chrl = [] # list of chromsome names.
        for n in xrange(self.nChromosomes):
            self.chrl.append(str(n))
        self.chrl.append("X")
        self.chrl.append("Y")

        # known gene list required pile-up code.
        #self.features.append(feature("known Gene", (os.path.join(self.path, "knownGene.txt")), "Gene")) # knownGne.txt
        self.features.append(feature("refseq", (os.path.join(self.path, "refseq.tsv")), "Gene")) # knownGne.txt

    def getListOfChromosomes(self):
        """
        return a list of the chromosome names
        (And thus corresponding keys)
        """
        return(self.chrl)

    def bindFeature(self, featureFile, featureType):
        """
        bind a 'feature file' to the genome
        """
        f = feature.feature(featureFile, featureType)
        self.features.append(f)

    def getSequenceA(self, chr, left, right):
        """
        extract the sequence using chr, left, right
        """
        pass

    def getSequenceB(self, location):
        """
        extract some sequence using chrX:left-right
        """
        pass

    def getFeatureA(self, chr, left, right, featureType):
        """
        get the feature spanning chr, left, right
        """
        pass

    def getFeatureB(self, location, featureType):
        """
        get the features spanning chr, left, right using chr:left-right notation
        """
        pass

    def getAllFeaturesInRange(self, location, feature=False):
        """
        get all of the draw features within the range location.
        overrides the feature version...
        """
        loc = utils.getLocation(location)

        # which feature to get?
        # for now just known gene:
        if feature:# I have a particule fature in mind:
            return([])
        else: # no idea, just return all;
            return(self.features[0].getAllFeaturesInRange(loc["chr"], loc["left"], loc["right"]))

if __name__ == "__main__":
    g = genome("mouse", "mm8")

