"""
genome pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
import glbase # import like this to get around namespace issues

def qcollide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return(False)
    if Aleft > Bright:
        return(False)

    if Aleft <= Bright and Aright >= Bright:
        return(True) # Bright point is within A, collision

    if Aright >= Bleft and Aleft <= Bleft:
        return(True) # Bleft point is within A, collision.

    if Bleft <= Aright and Bright >= Aright:
        return(True) # Aright point is within B, collision

    if Bright >= Aleft and Bleft <= Aleft:
        return(True) # Aleft point is within B, collision.

    #print "unhandled!"
    return(False)

class genome(glbase.genome):
    """
    chipFish override of vanilla glbase

    Adds descriptions for chipFish to interpret the interfaces from
    """

    # extra methods:
    def getAllDrawableFeaturesInRange(self, location):
        """
        (Extra)
        retrieve all features between location.

        (To move into glbase proper?)
        """
        print location
        ret = []
        if self.dataByChr.has_key(location["chr"]):
            for item in self.dataByChr[location["chr"]]:
                print location["left"], location["right"], item["loc"]["left"], item["loc"]["right"]
                if qcollide(location["left"], location["right"], item["loc"]["left"], item["loc"]["right"]):
                    print item
                    """
                    {'loc': location: contents: {'chr': '3', 'right': 153791001, 'left': 153781632},
                    'name': 'Asb17',
                    'n': 6544,
                    'array_systematic_name': 'scl22442.2.1_3-S',
                    'tss_loc': location: contents: {'chr': '3', 'right': 153781632, 'left': 153781632},
                    'refseq': 'NM_025758',
                    'entrez': 66772,
                    'strand': '+',
                    'description': 'ankyrin repeat and SOCS box-containing 17'}
                    """
                    # make a suitable draw object
                    item["type"] = "gene" # set the type flag for gDraw
                    ret.append(item)
        return(ret)

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

    # available for the gui on this class
    __gui__avail__ = {"load a list": glbase.genome.load,
        "append entry": glbase.genome.append,
        "collide lists": glbase.genome.collide,
        "overlap lists": glbase.genome.overlap,
        "map lists": glbase.genome.map,
        "remove duplicates": glbase.genome.removeDuplicates,
        "reverse": glbase.genome.reverse,
        "save list": glbase.genome.save,
        "sort list by key": glbase.genome.sort
        }
