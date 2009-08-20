"""

data.py

override for chipFish

contains misc data from glbase, stuff found in helpers, data, flags etc...

(c) 2009 oAxiom

"""
import sys, os

sys.path.append("..")
import opt

sys.path.append(opt.path.glbase_package)
import glbase # import like this to get around namespace issues

from glbase.flags import * # just pass through.
from glbase.helpers import load, fold2UpOrDown
