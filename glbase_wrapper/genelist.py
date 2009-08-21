"""
genelist pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
import glbase # import like this to get around namespace issues

class genelist(glbase.genelist):
    """
    chipFish override of vanilla glbase

    Adds descriptions for chipFish to interpret the interfaces from
    """
    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Vanilla genelist"

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {"load a list": glbase.genelist.load,
        "append entry": glbase.genelist.append,
        "collide lists": glbase.genelist.collide,
        "overlap lists": glbase.genelist.overlap,
        "map lists": glbase.genelist.map,
        "remove duplicates": glbase.genelist.removeDuplicates,
        "reverse": glbase.genelist.reverse,
        "save list": glbase.genelist.save,
        "sort list by key": glbase.genelist.sort
        }

    # I have to redefine the methods chipFish respects here:
    def annotate(self, **kargs):
        glbase.genelist.annotate(self, **kargs)

    # define gui interfaces
    annotate.tooltype = "Annotate a list of genes"
    annotate.__gui__ = {"required": {"list": "genelist", "key": "list key", "distance": "int"},
        "optional": {"resolution": "filename"}} # bind gui descriptor

    # next method:
