"""
constants, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for immutable constants

"""

import cairo

# genuine constants, unchangable:

# text constants
TXT_NORMAL = 1
TXT_BOLD = 2
TXT_ITALIC = 3

# convert my constants to cairo text constants.
txtToCairoA = {None: cairo.FONT_WEIGHT_NORMAL,
        TXT_NORMAL: cairo.FONT_WEIGHT_NORMAL,
        TXT_BOLD: cairo.FONT_WEIGHT_BOLD,
        TXT_ITALIC: cairo.FONT_WEIGHT_NORMAL
        }

txtToCairoB = {None: cairo.FONT_SLANT_NORMAL,
        TXT_NORMAL: cairo.FONT_SLANT_NORMAL,
        TXT_BOLD: cairo.FONT_SLANT_NORMAL,
        TXT_ITALIC: cairo.FONT_SLANT_ITALIC
        }

validChrNames = ["X", "Y", "M"]
