"""
chromosome, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

"""

import sys, os, csv

import error, db, utils

validFeatureTypes = [
    "Point", "Span", "Peak", "Expression", # generic types
    "Gene", "ncRNA", "histone" # specific architecture
    ]

class feature:
    def __init__(self, name, featureFile, featureType):
        """
        This is the implemetation of 'features'
        some sort of decorators of the genome.
        They are many types describing a variety of objects.
        """
        self.type = featureType
        self.file = featureFile
        self.name = name
        self.data = {} # arranged by Chr.

        self.importFeature(featureFile)

    def _exportFeatures(self):
        """
        Pickle the feature
        """
        pass

    def _importFeatures(self):
        """
        unPickle the feature
        """
        pass

    def drawFeature(self, location):
        """
        draw this feature across location to surface.
        interfaces with gDraw
        """
        pass

    def importFeature(self, filename):
        """
        import a feature from filename.
        this will interface with the db later, which will handle import/export.
        """
        oh = open(filename, "rU")
        reader = csv.reader(oh, dialect=csv.excel_tab)
        print self.data

        # currenlty a custom loader for knownGene. table.

        # [knownGene table layout]
        # 0name, 1chr, 2strand, 3start, 4end, 5cdssart, 6cdsend,
        # 7exonCount, 8exonStart, 9exonEnd, 10protienID, 11alignID

        for line in reader:
            if not line[0].count("#"):
            #{"type": "lncRNA", "strand": "-", "left": 10000, "right": 11000}
                d = {"type": "gene",
                    "strand": line[2],
                    "chr": line[1].strip("chr"),
                    "start": int(line[3]),
                    "end": int(line[4]),
                    "cdsStart": int(line[5]),
                    "cdsEnd": int(line[6]),
                    "exonCount": int(line[7]),
                    "exonStarts": [int(x) for x in line[8].strip(",").split(",")],
                    "exonEnds": [int(x) for x in line[9].strip(",").split(",")],
                    "proteinID": line[10],
                    #"alignID": line[11],
                    "name": line[11], # these two only in refseq list.
                    "locuslink": line[12]
                    } # this is how it will come out of the db.
                if not self.data.has_key(d["chr"]):
                    self.data[d["chr"]] = []
                self.data[d["chr"]].append(d)
        oh.close()
        return(True)

    def getAllFeaturesInRange(self, chr, left, right):
        """
        return all of the data within the range chr, left, right.
        """
        #print chr, left, right
        ret = []
        if self.data.has_key(chr):
            for item in self.data[chr]:
                #print left,right, item["left"], item["right"]
                if utils.collide(left, right, item["start"], item["end"]):
                    ret.append(item)
        return(ret)



if __name__ == "__main__":
    f = feature("Known Genes", "/home/andrew/chipfish/chipfish/data/mm/mm8/refseq.tsv", "Gene")
    #print f.data
    print f.data["3"][0:10]
    print f.getAllFeaturesInRange("3", 159700000, 159820000)
