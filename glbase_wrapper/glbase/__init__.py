"""

Initialise glbase, import all the libraries, set up the environment etc.

Requires:
* numpy
* matplotlib
* scipy
"""

import sys

#-----------------------------------------------------------------------
# Load all of the global configuration options.
try:
    import config
    from errors import LibraryNotFoundError
except:
    print "Error: Fatal - GLbase is not installed correctly, cannot find my own libraries"
    print "       Is the python 'sys.path' correct?"
    sys.exit() # no raise if I can't get errors, it's surely a fatal installation problem.

# ----------------------------------------------------------------------
# Test for availability of the core non-standard libs.
# These need to be available as the subsequent load/checking is weak/non-existant.

try:
    import numpy
    config.NUMPY_AVAIL = True
except Exception:
    raise LibraryNotFoundError, "Fatal - Numpy is not available or not installed"

try:
    import scipy
    config.SCIPY_AVAIL = True
except Exception:
    raise LibraryNotFoundError, "Fatal - Scipy is not available or not installed"

try:
    import matplotlib
    matplotlib.use("Agg") # cluster friendly!
    config.MATPLOTLIB_AVAIL = True
except Exception:
    raise LibraryNotFoundError, "Fatal - matplotlib not available or not installed"

try:
    import sklearn
    config.SKLEARN_AVAIL = True
except Exception:
    pass # pass silently as sklearn is currently optional.

# ----------------------------------------------------------------------
# Now import the rest of my libraries - assumes here they are available.
# If I can get config and errors then these are probably available too.

from helpers import * # naughty, brings in data.py, cPickle, math, utils, config, sys, os, extra naughty as I probably don't need all that lot anyway now.
from location import location
from genelist import genelist
from expression import expression
from genome import genome
from delayedlist import delayedlist
from glglob import glglob
from element import motif
from track import track
from flat_track import flat_track
from progress import progressbar
from pwm import pwm
from pwms import pwms
from ecrbase import ecrbase, tfbs_iter
from region import region
from realtime import realtime
from expression import expression
from logos import logo
from draw import draw
from format_container import fc
from fastq import fastq
import gldata
import utils
import format
import cmaps

from tools.seqToTrk import seqToTrk
from tools.wigstep_to_flattrack import wigstep_to_flat
from tools.gerp_to_flattrack import gerp_to_flat
from tools.bigwig_to_flattrack import bigwig_to_flat
from tools.rnaseq import rnaseqqc

def version():
    config.log.info("glbase - version: %s %s" % (config.version, config.DATE))
    config.log.info("The working directory is: '%s'" % (os.getcwd()))

# export all of the libraries, methods and helpers.
__all__ = ["genelist", "fastq", "expression", "genome", "track", "flat_track", "delayedlist", # primary objects
            "location", "pwm", "pwms", #accesory objects 
            "flags",  "format",
            "utils", "glload", "seqToTrk", "logo", 
            "glglob", "motif",  "wigstep_to_flat", "bigwig_to_flat",
            "rnaseqqc", "gldata",
            "gerp_to_flat", "draw", "fc",
            "progressbar", "ecrbase", "region", "realtime", "tfbs_iter",
            "strandSorter",
            "cmaps"] + dir(helpers)
            # in future I want to get rid of dir() and control what gets exported.
