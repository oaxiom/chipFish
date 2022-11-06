"""

Initialise the glbase_wrapper, import me.

override for chipFish

contains misc data from glbase, stuff found in helpers, data, flags etc...

(c) 2009-2017 oAxiom

"""
import sys, os

# get chipFish options
#sys.path.append("../")
#import opt

# get the actual glbase package location (probably ./glbase)
#sys.path.append(opt.path.glbase_package)
#import glbase # import like this to get around namespace issues

from . import glbase_stage
VERSION = glbase_stage.config.VERSION

# pass through directly from glbase
from .glbase_stage.flags import * # just pass through.
from .glbase_stage.location import location
from .glbase_stage.helpers import * # pass through
from .glbase_stage.delayedlist import delayedlist
from .glbase_stage.tools.seqToTrk import seqToTrk
from .glbase_stage import format

# get the wrapped versions.
from .genelist import genelist
from .genome import genome
from .genome_sql import genome_sql
from .flat_track import flat_track

# make these functions available in the package
__all__ = ["genelist", "location", "genome", "delayedlist", "VERSION",
    "format",
    ] + dir(glbase_stage.flags) + dir(glbase_stage.helpers)
