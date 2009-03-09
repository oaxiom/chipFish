"""
opt, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for mutable constants and other data

. These should be loadable/savable and changable.
"""

# set of options

import sys, os

# generic options
app_path = sys.path[0]

# graphic options
class graphics:
    """
    general graphics options.
    """
    gene_height = 10        # width to draw the genes
    lncrna_height = 8
    microRNA_height = 6

    screen_colour = (1,1,1,1)
    gene_colour = (0,0,0)
    lncrna_colour = (0.1, 0.8, 0.1)
    microRNA_colour = (0.6, 0.1, 0)

class ruler:
    """
    options specific to the ruler.
    """
    ruler_font = "Arial"
    ruler_height_px = 10 # height ofruler in pixels.
    ruler_text_height = 20 # distance from top of screen of the ruler text.
    ruler_colour = (0,0,0)
    ruler_line_width = 1.5

if __name__ == "__main__":
    # error checking
    #print dir(graphics.gene_colour)
    print graphics.gene_colour

