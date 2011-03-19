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
    _available_draw_types = ("bar", "spot")

    def __init__(self, *args, **kargs):
        glbase.peaklist.__init__(self, *args, **kargs)
        
        self.meta_data = {"name": self.name}

    def __getitem__(self, key):
        """
        Emulate a dict
        """
        if key == "info": # catch this special key
            for k in self.meta_data:
                print "%s\t:\t%s" % (k, self.meta_data[k])
        else:
            assert key in self.meta_data, "'%s' not found in this track" % key
            return(self.meta_data[key])

    # new methods:
    def get_data(self, type, loc, strand=None, resolution=1, **kargs):
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
        if type == "bar":
            # make a single array
            ret = zeros(int( (loc["right"]-loc["left"]+resolution)/resolution ))

            if loc["chr"] in self.dataByChr:
                for item in self.dataByChr[loc["chr"]]:
                    if qcollide(loc["left"], loc["right"], item["loc"]["left"], item["loc"]["right"]):
                        # make a suitable draw object
                        for rloc in xrange(item["loc"]["left"], item["loc"]["right"], int(resolution)):
                            array_relative_location = int((rloc - loc["left"]) / resolution) # convert relative to the array

                            if array_relative_location >= 0 and array_relative_location < len(ret): # within array
                                ret[array_relative_location] += 1
            return(ret)
        elif type == "spot":
            ret = []
            if loc["chr"] in self.dataByChr:
                for item in self.dataByChr[loc["chr"]]:
                    if qcollide(loc["left"], loc["right"], item["loc"]["left"], item["loc"]["right"]):
                        ret.append(item["loc"])
            return(ret)
        return(None)

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        "draw type": _available_draw_types
        }


