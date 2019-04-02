"""
genome pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009-2019 oAxiom


"""

import sys, os

from .utils import qcollide

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
from . import glbase_stage # import like this to get around namespace issues

class genome_sql(glbase_stage.genome_sql):
    """
    chipFish override of vanilla glbase

    Adds descriptions for chipFish to interpret the interfaces from
    """
    # descriptions, text for localisation
    name = 'None' # For compatability with vanilla genomes
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Genome"
    _default_draw_type = "repeats"
    _available_draw_types = ("genome", "genome_stack") # genome_stack not implemented.

    # extra methods:
    def get_data(self, location):
       return(self.getFeatures(location))

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
        glbase.genome.annotate(self, **kargs)

    # define gui interfaces
    annotate.tooltype = "Annotate a list of genes"
    annotate.__gui__ = {"required": {"list": "genelist", "key": "list key", "distance": "int"},
        "optional": {"resolution": "filename"}} # bind gui descriptor

    # next method:
