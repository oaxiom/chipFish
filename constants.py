"""
constants, part of chipFish

(c) 2008-2012 oAxiom

Not for distribution.

Class container for immutable constants

"""

import cairo

# convert my constants to cairo text constants.
txtToCairoA = {None: cairo.FONT_WEIGHT_NORMAL,
        "normal": cairo.FONT_WEIGHT_NORMAL,
        "bold": cairo.FONT_WEIGHT_BOLD,
        "bolditalic": cairo.FONT_WEIGHT_BOLD,
        "italic": cairo.FONT_WEIGHT_NORMAL
        }

txtToCairoB = {None: cairo.FONT_SLANT_NORMAL,
        "normal": cairo.FONT_SLANT_NORMAL,
        "bold": cairo.FONT_SLANT_NORMAL,
        "bolditalic": cairo.FONT_SLANT_ITALIC,
        "italic": cairo.FONT_SLANT_ITALIC
        }

validChrNames = frozenset(["X", "Y", "M"])
