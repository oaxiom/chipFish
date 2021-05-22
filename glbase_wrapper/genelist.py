"""
genelist pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os, numpy

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
from . import glbase_stage # import like this to get around namespace issues

class genelist(glbase_stage.genelist):
    """
    chipFish override of vanilla glbase

    Adds descriptions for chipFish to interpret the interfaces from
    """
    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Vanilla genelist"
    _default_draw_type = "spot"
    _available_draw_types = ("bar", "spot")

    # I have to redefine the methods chipFish respects here:
    def annotate(self, **kargs):
        glbase_stage.genelist.annotate(self, **kargs)

    # define gui interfaces
    annotate.tooltype = "Annotate a list of genes"
    annotate.__gui__ = {"required": {"list": "genelist", "key": "list key", "distance": "int"},
        "optional": {"resolution": "filename"}} # bind gui descriptor

    def __getitem__(self, key):
        """
        Emulate a dict

        I need to override this, to only return metadata.
        glbase genome objects return self.qkeyfind[key] by default.
        """
        if key == "info": # catch this special key
            for k in self.meta_data:
                print(("%s\t:\t%s" % (k, self.meta_data[k])))
        elif key == "name":
            return(self.name)
            #assert key in self.meta_data, "'%s' not found in this track" % key
            #return(self.meta_data[key])

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
        ret = []
        if loc["chr"] in self.dataByChr:
            for item in self.dataByChr[loc["chr"]]:
                if loc.qcollide(item["loc"]):
                    ret.append(item["loc"])
        return(ret)

