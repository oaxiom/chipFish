"""
opt, part of chipFish

(c) 2008-2015 oAxiom

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
    huge_move = 100
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
    repeat_height = 10 # Height of the repeats

    screen_colour = (1,1,1,1)
    gene_colour = (0,0,0) #Black
    lncrna_colour = (0.1, 0.8, 0.1)
    microRNA_colour = (0.6, 0.1, 0)
    repeat_cols = {"LINE": (0, 0, 0.7), "SINE": (0,0.7,0), "LTR": (0.7,0,0), "DNA": (0.5, 0.5, 0.5)}
        
    # Font related defaults
    font = "Helvetica"
    default_font_size = 11 # Fall back font size
    repeat_label_font_size = 12
    repeat_label_font_style = 'normal'

    right_border_width = 0 # in pixels, size of rightmost click border.

    # gene arrow
    arrow_width_px = 4 # For genes
    arrow_height_px = 4
    repeat_arrow_height_px = 6 # Fore repeats

class draw:
    """
    drawing options
    """
    double_lines_for_genome = False
    single_midline_in_introns = True # draws a line through the gene, enclosing the gene
    chevrons_inside_introns = False # not implemented draws chevrons inside the gene indicating direction
    braces_between_exons = False # not implemented draws braces across the exons
    scale_bar = True # The scale bar at the top right of the image just below the ruler
    scale_bar_fontsize = 16 # scale bar font size. 
    view_port_width = 1000 # The width of the genome view in pixels
    genomic_location = True # draw the chrN:nnnnnn-nnnnnn label just above the first track
    genomic_location_font_size = 28 # font size for the genomic location label

class gene:
    font_size = 26
    font_style = "normal"

class ruler:
    """
    options specific to the genome ruler.
    """
    font_size = 18
    height_px = 10 # height of ruler in pixels.
    text_height = 10 # distance from top of screen of the ruler text.
    colour = (0,0,0)
    line_width = 1
    draw = False

class track:
    """
    options for the tracks
    """
    height_px = {"graph": 200, # These are the pixel heights of the tracks. Can be any reasonable value. These values are asthetically pleasing
        "graph_split_strand": 200, # same as graph
        "kde_graph": 200, # ditto
        "bar": 20,
        "spot": 20,
        "genome": 70,
        "repeats": 200} # height's of the tracks in pixels.
    bar_height = 12 # Height of the bars for bar tracks    
    spot_pixel_radius = 7 # Size of the spot circle
    spot_default_colour = (0.1, 0.1, 0.1) # The colour of the spot for spot tracks
    bar_default_colour = (0.1, 0.1, 0.1) # The colour of the spot for spot tracks
    spot_filled = True # fill the spot circle?
    spot_shape = "circle" # supported = circle, triangle
    font_scale_size = 24 # Size of the font on the tracks
    filled = True
    background = False # Draw a grey background defining the span of the track
    draw_names = True
    draw_scales = True
    scale_bar_font_size = 52
    label_fontsize = 22
    min_scale = 20 # The minimum x axis value for tracks. 
    lock_scales = False # Lock all scales on the Tracks together. Make this per trk file settable?