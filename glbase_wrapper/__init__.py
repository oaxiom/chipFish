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
from glbase.config import VERSION
sys.path.append("glbase/tools/") # seqToTrk is not exported by glbase in rev 131.
#from seqToTrk import seqToTrk # throws an error for some reason?

# get the wrapped versions.
from genelist import genelist
from genome import genome
from peaklist import peaklist
from track import track

# make these functions available in the package
__all__ = ["genelist", "location", "genome", "delayedlist", "VERSION",
    "track", "peaklist"#"seqToTrk"
    ] + dir(glbase.flags) + dir(glbase.helpers) # in future I want to get rid of dir() and control what gets exported.
