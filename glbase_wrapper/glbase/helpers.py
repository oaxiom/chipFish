"""

helpers.py

a selection of functions to help various parts of the code.

Used in places like the CSV format language and in selecting criteria, etc...

"""

import cPickle, math, sys, os
import utils, config

from data import *
from location import location

# ----------------------------------------------------------------------
# Some helper functions

def glload(filename):
    """
    load a pickle.
    """
    assert os.path.exists(os.path.realpath(filename)), "File '%s' not found" % filename

    oh = open(os.path.realpath(filename), "rb")
    newl = cPickle.load(oh)
    oh.close()
    config.log.info("Loaded '%s' binary file with '%s' items" % (filename, len(newl)))
    return(newl)

def change_drawing_mode(self, mode):
    """
    Change the output driver

    send me "png", "eps"|"ps" for your output needs
    """
    recognised_modes = ["png", "ps", "eps"]
    self.output_mode = "png"

    assert mode in recognised_modes, "Draw output mode '%s' is not recognised" % mode

    if mode == "png":
        config.DEFAULT_DRAWER = "png"
    elif mode == "eps" or mode == "ps":
        config.DEFAULT_DRAWER = "eps"
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

# various other helpers for normalisation etc..
