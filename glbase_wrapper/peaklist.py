"""
peaklist pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os

from numpy import zeros

from utils import qcollide

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
import glbase # import like this to get around namespace issues

class peaklist(glbase.peaklist):
    """
    chipFish override of vanilla glbase track.py

    In chipFish I want to store some information about the current track
    and
    """
    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Vanilla peaklist"
    _default_draw_type = "bar"

    # new methods:
    def get_array(self, loc, strand=None, resolution=1, **kargs):
        """
        **Purpose**
            mimic the entry to track.

        **Arguments**
            loc
                location span to get the array for.

            strand (Not supported)

            resolution
                the bp resolution to use for the array.

        **Results**
            returns a Numpy array.
        """
        # make a single array
        a = zeros(int( (loc["right"]-loc["left"]+resolution)/resolution ))

        ret = []
        if loc["chr"] in self.dataByChr:
            for item in self.dataByChr[loc["chr"]]:
                if qcollide(loc["left"], loc["right"], item["loc"]["left"], item["loc"]["right"]):
                    # make a suitable draw object
                    item["type"] = "bar" # set the type flag for gDraw

                    for rloc in xrange(item["loc"]["left"], item["loc"]["right"], int(resolution)):
                        array_relative_location = int((rloc - loc["left"]) / resolution) # convert relative to the array

                        if array_relative_location >= 0 and array_relative_location < len(a): # within array
                            a[array_relative_location] += 1
        return(a)

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        "draw type": {"type": "selection", "options": ("bar", "spots")}
        }


