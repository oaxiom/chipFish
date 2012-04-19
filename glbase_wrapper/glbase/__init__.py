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
    config.MATPLOTLIB_AVAIL = True
except Exception:
    raise LibraryNotFoundError, "Fatal - matplotlib not available or not installed"

# ----------------------------------------------------------------------
# Now import the rest of my libraries - assumes here they are available.
# If I can get config and errors then these are probably available too.

from flags import *
from helpers import *
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
from draw import draw # draw is available?
from format_container import fc
import utils
import format

from tools.seqToTrk import seqToTrk
from tools.wigstep_to_flattrack import wigstep_to_flat
from tools.gerp_to_flattrack import gerp_to_flat

config.log.info("glbase - version: %s %s" % (config.version, config.DATE))
config.log.info("The working directory is: '%s'" % (sys.path[0]))

# export all of the libraries, methods and helpers.
__all__ = ["genelist", "flags", "expression", "genome", "format",
            "utils", "glload", "seqToTrk", "logo", "delayedlist",
            "glglob", "motif", "track", "flat_track", "wigstep_to_flat", 
            "gerp_to_flat", "draw", "fc",
            "progressbar", "ecrbase", "region", "realtime", "tfbs_iter",
            "location", "pwm", "pwms", "expression", "glload", "strandSorter"] + dir(helpers)
            # in future I want to get rid of dir() and control what gets exported.
