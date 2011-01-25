"""

Initialise glbase, import all the libraries, set up the environment etc.

Requires:
* numpy
* matplotlib
* scipy
"""

import sys, os

#-----------------------------------------------------------------------
# Load all of the global configuration options.
try:
    import config
    from errors import LibraryNotFoundError
except:
    print "Error: Fatal - GLbase is not installed correctly, cannot find my own libraries"
    print "       Is the python 'sys.path' correct?"
    sys.exit() # no raise if I can't get errors.

# ----------------------------------------------------------------------
# Test for availability of the core non-standard libs.
# These need to be available as the subsequent load/checking is weak/non-existant.

try:
    import numpy
    config.NUMPY_AVAIL = True
except:
    raise LibraryNotFoundError, "Fatal - Numpy is not available or not installed."

try:
    import scipy
    config.SCIPY_AVAIL = True
except:
    raise LibraryNotFoundError, "Fatal - Scipy is not available or not installed."

try:
    import matplotlib
    config.MATPLOTLIB_AVAIL = True
except:
    raise LibraryNotFoundError, "Fatal - matplotlib not available or not installed."

# ----------------------------------------------------------------------
# Now import the rest of my libraries - assumes here they are available.

from flags import *
from helpers import *
from location import location
from genelist import genelist
from microarray import microarray
from genome import genome
from taglist import taglist
from delayedlist import delayedlist
from peaklist import peaklist
from glglob import glglob
from element import motif
from track import track
from progress import progressbar
from pwm import pwm
from pwms import pwms
from ecrbase import ecrbase, tfbs_iter
from region import region
from realtime import realtime
import draw # draw is available?
import utils

from tools.seqToTrk import seqToTrk

config.log.info("glbase - version: %s, %s" % (config.VERSION, config.DATE))

# export all of the libraries, methods and helpers.
__all__ = ["genelist", "flags", "microarray", "genome",
            "taglist", "utils", "glload", "seqToTrk",
            "peaklist", "glglob", "motif", "track",
            "progressbar", "ecrbase", "region", "realtime", "tfbs_iter",
            "location", "pwm", "pwms"] + dir() # in future I want to get rid of dir() and control what gets exported.

