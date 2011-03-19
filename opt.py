"""
opt, part of chipFish

(c) 2008-2011 oAxiom

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
    web = True
    app_path = sys.path[0]
    v = ["0.3"]
    if web:
        v.append("web")
    if debug:
        v.append("debug")
        
    version = " ".join(v)

__VERSION__ = generic.version

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
    huge_move = 90
    big_move = 30 # in percent
    small_move = 10 # in percent

class debug:
    """
    containers for debugging options.
    """
    profile = False
    draw_collision_boxes = False

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

    # Font related defaults
    font = "sans-serif"
    default_font_size = 12 # Generally refers here to the track Text sizes

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
    font = "sans-serif"
    font_size = 10
    height_px = 10 # height of ruler in pixels.
    text_height = 10 # distance from top of screen of the ruler text.
    colour = (0,0,0)
    line_width = 1

class track:
    """
    options for the tracks
    """
    height_px = {"graph": 100,
        "graph_split_strand": 100, # same as graph
        "kde_graph": 100, # ditto
        "bar": 30,
        "spot": 30,
        "genome": 50} # height's of the tracks in pixels.
    spot_pixel_radius = 4
    spot_default_colour = (0.8, 0.1, 0.1)
    spot_filled = True # fill the spot circle
    spot_shape = "circle" # supported = circle, triangle
    font_scale_size = 8