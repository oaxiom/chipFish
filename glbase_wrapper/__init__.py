"""

Initialise the glbase_wrapper, import me.

override for chipFish

contains misc data from glbase, stuff found in helpers, data, flags etc...

(c) 2009 oAxiom

"""
import sys, os

# get chipFish options
sys.path.append("..")
import opt

# get the actual glbase package location (probably ./glbase)
sys.path.append(opt.path.glbase_package)
import glbase # import like this to get around namespace issues

# pass through directly from glbase
from glbase.flags import * # just pass through.
from glbase.location import location
from glbase.helpers import * # pass through
from glbase.delayedlist import delayedlist

# get the wrapped versions.
from genelist import genelist
from genome import genome


# make these functions available in the package
__all__ = ["genelist", "location", "genome", "delayedlist",
    ] + dir(glbase.flags) + dir(glbase.helpers) # in future I want to get rid of dir() and control what gets exported.
