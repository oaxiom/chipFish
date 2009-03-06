"""
genome, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

"""

import feature

import chromosome as chr

spcode = {
    "Mus Musculus": "MM",
    "mouse": "MM",
    "Homo Sapiens": "HS",
    "human": "HS"
    }

class genome:
    def __init__(self, species):
        """
        Initialise the genome,
        species is a simple species string that corresponds to one of the species codes.
        """
        self.chrl = []
        self.species = species
        self.features = []
        if species.lower() in spcode:
            self.code = spcode[species.lower()]
        else:
            self.code = "UU"

    def addChromosome(self, number):
        c = chr()
        c.setNumber(number)
        self.chrl.append(c)

    def bindFeature(self, featureFile, featureType):
        """
        bind a 'feature file' to the genome
        """
        f = feature.feature(featureFile, featureType)
        self.features.append(f)

    def getSequenceA(self):
        """
        extract the sequence using chr, left, right
        """
        pass

    def getSequenceB(self):
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
