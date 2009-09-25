"""
opt, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for mutable constants and other data

. These should be loadable/savable and changable.
. save as a editable text file?
"""

import sys, os

# generic options
class generic:
    app_path = sys.path[0]

class path:
    glbase_location = os.path.realpath("./glbase_wrapper/glbase")
    glbase_package = os.path.realpath("./glbase_wrapper/")
    glbase_wrapper = os.path.realpath(".")

class debug:
    """
    containers for debugging options.
    """
    debug_version = False
    if debug_version:
        draw_collision_boxes = True
    else:
        draw_collision_boxes = False

class graphics:
    """
    general graphics options.
    """
    gene_height = 8        # width to draw the genes
    cds_height = 12
    lncrna_height = 7
    microRNA_height = 6

    screen_colour = (1,1,1,1)
    gene_colour = (0.89,0.654,0.165) #E6A729
    lncrna_colour = (0.1, 0.8, 0.1)
    microRNA_colour = (0.6, 0.1, 0)

    font = "Arial"

class ruler:
    """
    options specific to the genome ruler.
    """
    font = "Arial"
    font_size = 6
    height_px = 10 # height of ruler in pixels.
    text_height = 10 # distance from top of screen of the ruler text.
    colour = (0,0,0)
    line_width = 1

def saveOptions():
    pass

def loadOptions():
    pass


if __name__ == "__main__":
    # error checking
    #print dir(graphics.gene_colour)
    # export the options as a csv.
    pass

