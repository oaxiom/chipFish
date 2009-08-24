"""

helpers.py

a selection of functions to help various parts of the code.

Used in places like the CSV format language and in selecting criteria, etc...

"""

import cPickle, math, sys, os
import utils, config

from location import location

positive_strand_labels = frozenset(["+", "1", "f", "F"])
negative_strand_labels = frozenset(["-", "0", "r", "R"])

# ----------------------------------------------------------------------
# Some helper functions

def load(filename):
    """
    load a pickle.
    """
    oh = open(filename, "rb")
    newl = cPickle.load(oh)
    oh.close()
    return(newl)

def change_drawing_mode(self, mode):
    """
    Change the output driver

    send me "png", "eps"|"ps" for your output needs
    """
    recognised_modes = ["png", "ps", "eps"]
    self.output_mode = "png"

    assert mode in recognised_modes, "Draw output mode %s is not recognised" % mode

    if mode == "png":
        config.DEFAULT_DRAWER = "png"
    elif mode == "eps" or mode == "ps":
        config.DEFAULT_DRAWER = "eps"
    else:
        print "Error: Draw Output mode: %s not recognised" % mode
        config.DEFAULT_DRAWER = "png"
        return(False)
    return(True)

# ----------------------------------------------------------------------
# Criteria functions.

def fold2UpOrDown(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data > 2:
                    return(True)
                elif normed_data < 0.5:
                    return(True)
                else:
                    return(False)
        return(False)

def fold2Down(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data < 0.5:
                    return(True)
                else:
                    return(False)
        return(False)

def fold2Up(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data > 2:
                    return(True)
                else:
                    return(False)
        return(False)

def XDown(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the fold2up/down data
    """
    X = kargs["X"]
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data < X:
                    return(True)
                else:
                    return(False)
        return(False)

def XUp(data, names, normed = None, **kargs):
    """
    good for normalised illumina data, returns the Xdown data
    """
    X = kargs["X"]
    if normed:
        norm_value = data[normed]
        for c in data:
            if c != normed:
                normed_data = (data[c] / data[normed])
                # this is greedy - only 1 condition needs to fulfill the criteria.
                if normed_data > X:
                    return(True)
                else:
                    return(False)
        return(False)

# For formatting the CSV loading.

def strandSorter(chr, left, right, strand):
    """
    A helper proc to extract the tss from a list of coords.
    """
    if strand in positive_strand_labels:
        return(location(chr=chr, left=left, right=left))
    elif strand in negative_strand_labels:
        return(location(chr=chr, left=right, right=right))
    return(None)

def barSplitter(value):
    return(value.split("|"))

# various other helpers for normalisation etc..
