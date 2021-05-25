"""
genome pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os

from .utils import qcollide

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
from . import glbase_stage # import like this to get around namespace issues

class genome(glbase_stage.genome):
    """
    chipFish override of vanilla glbase

    Adds descriptions for chipFish to interpret the interfaces from
    """
    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Genome"
    _default_draw_type = "genome"
    _available_draw_types = ("genome", "genome_stack") # genome_stack not implemented.


    # extra methods:
    def get_data(self, type, location):
        if type == "genome":
            return self.getFeatures(location)

    def getAllDrawableFeaturesInRange(self, location):
        """
        Historical place-holder, for deletion
        """
        return(self.get_data("genome", location))

    def find(self, text, guage_object=None, case_sensitive=False):
        """
        **Purpose**
            find any items containing 'text' in the genelist

        **Arguments**
            text (string, required)
                the text string to search for

            guage_object
                a progress bar of some description with the
                method SetRange().
                If None, the progressbar is ignored

            case_sensitive (Optional, default=False)
                treat the search as case senesitive (or not)

        **Returns**
            A list of dicts containing the search string
            modify a progress bar (NotImplemented)

        """
        ret = []

        if guage_object:
            guage_object.SetRange(len(self))

        for index, item in enumerate(self.linearData):
            for key in item:
                # can only search strings:
                tstr = str(item[key]) if case_sensitive else str(item[key]).lower()
                if text in tstr:
                    ret.append(item)
                    break # only append once
            if guage_object:
                guage_object.SetValue(index)

        return(ret)

    def __getitem__(self, key):
        """
        I need to override this, to only return metadata.
        glbase genome objects return self.qkewyfind[key] by default.
        """
        if key == "name":
            return(self.name)

    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Vanilla genelist"

    # gui stuff.

    # I have to redefine the methods chipFish respects here:
    def annotate(self, **kargs):
        glbase_stage.genome.annotate(self, **kargs)

    # define gui interfaces
    annotate.tooltype = "Annotate a list of genes"
    annotate.__gui__ = {"required": {"list": "genelist", "key": "list key", "distance": "int"},
        "optional": {"resolution": "filename"}} # bind gui descriptor

    # next method:
