"""

glglobs are almagamated lists, useful for drawing comparisons between lists.

glglobs also inherit from genelist, but many of the methods are broken.

renamed to glglob as it clashes with a matplotlib module name

"""

import sys, os, csv, string, math

from numpy import array
from array import array as qarray
from copy import deepcopy

import config
from flags import *
from genelist import genelist
from draw import draw
from errors import AssertionError
if config.MATPLOTLIB_AVAIL:
    import matplotlib.pyplot as plot

class glglob(genelist): # cannot be a genelist, as it has no keys...
    def __init__(self, *args, **kargs):
        genelist.__init__(self)

        # args should be a list of lists.
        # we then store them in the linearData set
        for a in args:
            for b in args:
                assert a == b, "Lists must be equivalent"

        self.linearData = args
        self._optimiseData()

    def __repr__(self):
        return("glbase.glglob")

    def _optimiseData(self): # no keys, so would die.
        pass

    def loadCSV(self):
        print "Error: glglobs cannot be sensibly represented as a CSV file, use load() to save binary copies"
        return(False)

    def saveCSV(self):
        print "Error: glglobs cannot be sensibly represented as a CSV file, use save() to save binary copies"
        return(False)

    def compare(self, key=None, filename="", **kargs):
        assert filename, "Filename must be a valid string"
        assert self.linearData[0].linearData[0].has_key(key), "Key must be valid" # just check one list

        print "Info: This may take a while, based upon the sizes of the lists. We have to intersect them all by %s" % key

        matrix = [[0 for l in xrange(self)] for l in xrange(self)] # 2D matrix.

        for ia, la in enumerate(self.linearData):
            for ib, lb in enumerate(self.linearData):
                if ia == ib:
                    matrix[ia][ib] = 0.0
                else:
                    if key == "loc":
                        matrix[ia][ib] = len(la.collide(gene_list=lb, loc_key=key, delta=200))
                    else:
                        matrix[ia][ib] = la.map(gene_list=lb, key=key)
        # normalise the matrix

        # draw a heatmap.

        return(True)


