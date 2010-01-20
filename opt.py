"""
opt, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for mutable constants and other data

. These should be loadable/savable and changable.
. save as a editable text file (No, this is the editable file)!
. But perhaps an editor is required? - If I get funding!
. This is actually more accurately consts and options together?
"""

import sys, os

# generic options
class generic:
    """
    some very generic options and constants.
    """
    debug = True
    app_path = sys.path[0]
    VERSION = "0.1"

class path:
    """
    some common path locations.
    """
    raw_glbase_location = os.path.realpath("./glbase_wrapper/glbase")
    glbase_package = os.path.realpath("./glbase_wrapper/")
    glbase_wrapper = os.path.realpath(".")

class interface:
    """
    user interface options
    """
    big_move = 30 # in percent
    small_move = 10 # in percent

class debug:
    """
    containers for debugging options.
    """
    profile = False
    draw_collision_boxes = True

class graphics:
    """
    graphics options.
    """
    gene_height = 8        # width to draw the genes
    cds_height = 12
    lncrna_height = 7
    microRNA_height = 6

    screen_colour = (1,1,1,1)
    gene_colour = (0.89,0.654,0.165) #E6A729
    gene_colour = (0,0,0) #Black
    lncrna_colour = (0.1, 0.8, 0.1)
    microRNA_colour = (0.6, 0.1, 0)

    font = "Arial"

    right_border_width = 20 # in pixels, size of rightmost click border.

    # gene arrow
    arrow_width_px = 4
    arrow_height_px = 4

class draw:
    """
    drawing options
    """
    double_lines_for_genome = False
    single_midline_in_introns = True # draws a line through the gene, enclosing the gene
    chevrons_inside_introns = False # not implemented draws chevrons inside the gene indicating direction
    braces_between_exons = False # not implemented draws braces across the exons
    scale_bar = True

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

class track:
    """
    options for the tracks
    """
    height_px = {"graph": 150,
        "graph_split_strand": 150, # same as graph
        "bar": 30,
        "spot": 30} # height's of the tracks in pixels.
    genome_base_offset = 60 # offset from the genome track to the first track.
    spot_pixel_radius = 4
    spot_default_colour = (0.8, 0.1, 0.1)
    spot_filled = True # fill the spot circle
    spot_shape = "circle" # supported = circle, triangle

def saveOptions():
    pass

def loadOptions():
    pass
