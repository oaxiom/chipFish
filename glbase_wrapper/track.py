"""
track pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
import glbase # import like this to get around namespace issues

class track(glbase.track):
    """
    chipFish override of vanilla glbase track.py

    In chipFish I want to store some information about the current track
    and
    """
    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Vanilla track"
    _default_draw_type = "graph"
    _available_draw_types = ("graph", "bar")

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        "draw type": _available_draw_types
        }


