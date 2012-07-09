"""

Initialise the glbase_wrapper, import me.

override for chipFish

contains misc data from glbase, stuff found in helpers, data, flags etc...

(c) 2009-2011 oAxiom

"""
import sys, os

# get chipFish options
#sys.path.append("../")
#import opt

# get the actual glbase package location (probably ./glbase)
#sys.path.append(opt.path.glbase_package)
#import glbase # import like this to get around namespace issues

import glbase
VERSION = glbase.config.VERSION

# pass through directly from glbase
from glbase.flags import * # just pass through.
from glbase.location import location
from glbase.helpers import * # pass through
from glbase.delayedlist import delayedlist
from glbase.tools.seqToTrk import seqToTrk 
import glbase.format as format

# get the wrapped versions.
from genelist import genelist
from genome import genome
from track import track
from flat_track import flat_track

# make these functions available in the package
__all__ = ["genelist", "location", "genome", "delayedlist", "VERSION",
    "track", "format",
    ] + dir(glbase.flags) + dir(glbase.helpers) 