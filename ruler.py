"""

ruler.py

abstraction of a generic ruler class.

"""

from __future__ import division
from error import AssertionError

import opt
from constants import *

class ruler:
    def __init__(self, min, max, pixel_span, units=None, useSI_label=True):
        """
        **Purpose**
            an abstract class to represent rulers.

            Only integer rulers are supported.

        **Arguments**
            min
                the minimum value of the ruler

            max
                the maximum value of the ruler.

            pixel_span (tuple)
                the size in pixels of the ruler.
                Currently non-zero starting pixel_spans are ignored

            useSI_label (default = True)
                append k for 1000, M for 1,000,000, G for 1,000,000,000 etc...
                to the labels.

        **Results**
            A valid ruler object
        """
        self.set(int(min), int(max), pixel_span)
        self.use_si = useSI_label
        self.units = units

    def set(self, min, max, pixel_span=None):
        """
        **Purpse**

        **Arguments**
            min
                update the minimum value

            max
                update the maximum value
                
            display_pixel_size (tuple)
                change the size in pixels of the viewport.
                leave as None to leave unchanged.
        """
        print pixel_span
        
        assert min < max, "min must be less than max for ruler"
        self.min = min
        self.max = max
        self.wid = self.max - self.min
        if pixel_span:
            self.display_px_w = pixel_span[1] - pixel_span[0] # Implemented, but probably not correct
            self.display_px_l = pixel_span[0] 
            self.display_px_r = pixel_span[1] # Currently ignored
        self.__units_in_px = ((self.min - self.max) / self.display_px_w) # number of units per pixel

    def get_ruler_data(self, percent=20, minor_percent=2):
        """
        **Purpose**

        **Arguments**

            percent (Optional)
                the percent to use to break up the ruler.
                (defaults to 20%)
                
            minor_percent (Optional)
                the percentage for the mintor ticks.

        **Returns**
        sends back a list containing tuple pairs in the form:
        [(location_px, "label"), (..) .. (..)]
        """
        fp = percent / 100.0

        ppx = round(self.display_px_w * fp) # work out how many pixels 20% is.
        p = round(self.wid * fp) # work out how many pixels 20% is.

        #print self.min, self.max, self.display_px_w
        #print "20%:", p, ppx, fp

        best = 0
        minim = 1e13

        # find the closest unit to p:
        scales_and_labels = {
            0.001: "m",
            0.01: "m",
            0.1: "m",
            1: "",
            10: "",
            100: "",
            1e3: "k",
            10e3: "k",
            100e3: "k",
            1e6: "M",
            10e6: "M",
            100e6: "M",
            1e9: "G",
            10e9: "G",
            100e9: "G", 
            1e12: "P" # As big as we go, should be okay for most purposes. a maxint is 1e16 anyway, if you are that big
            }         # Then this is the wrong module for you!
            # a 64-bit maxint:
            #9223372036854775807 = 1e16?
            
        best_div_lookback = {"m": 0.001, "": 1, "k": 1e3, "M": 1e6, "G": 1e9, "P": 1e12}
        # Look back for the number to divide by to maintain units.

        for i, v in enumerate(scales_and_labels):
            if minim > abs(p - v):
                minim = abs(p - v)
                best = v

        print "best:", best, scales_and_labels[best]
        
        major = []
        mi = int((self.min / int(best))) * int(best) # get the closest min
        ma = int((self.max / int(best))) * int(best) # and max
        for v in xrange(mi, ma, int(best)): 
            px_location = (self.display_px_w * ((v - self.min) / self.wid))
            if self.units:
                lab = "%s %s%s" % (int(v/best_div_lookback[scales_and_labels[best]]), scales_and_labels[best], self.units)
            else:
                lab = "%s %s" % (int(v/best), scales_and_labels[best])
            major.append( (int(px_location), lab.strip()) )

        return({"major": major, "minor": None})

    def draw(self, cairo_context, draw_top_left_tuple, label=None):
        """
        **Purpose**
            draw the ruler on a valid cairo_context device
            This is pretty internal to chipfish...

        **Arguments**
            cairo_context

            draw_top_left_tuple

        **Returns**
            a fourple containing the collision box coordinates
        """
        #x,y,w,h = self.__getTextExtents("Chromosome %s" % str(self.chromosome))
        #self.__drawText(5, opt.ruler.height_px + 22, opt.ruler.font, "Chromosome %s:%s-%s" % (str(self.chromosome), self.lbp, self.rbp), size=13)

        #self.__setPenColour(opt.ruler.colour)
        # work out a good scale representation
        # wether to draw at 100, 1000, 10000, 100000, 1000000 ...

        # current view delta = self.delta
        data = self.get_ruler_data()
        for item in data["major"]:
            print item
            cairo_context.set_line_width(opt.ruler.line_width)
            cairo_context.move_to(item[0], 0)
            cairo_context.line_to(item[0], opt.ruler.height_px)
            cairo_context.stroke()
            # draw a label
            
            cairo_context.set_source_rgb(0, 0, 0)
            cairo_context.select_font_face("sans-serif", cairo.FONT_WEIGHT_NORMAL, cairo.FONT_SLANT_NORMAL)
            cairo_context.set_font_size(8)
            cairo_context.move_to(item[0]+2, opt.ruler.height_px)
            cairo_context.show_text(str(item[1]))
        return(None)

if __name__ == "__main__":
    #min, max, display_pixel_size, useSI_label=True
    r = ruler(5, 10, (0, 100), "bp", True)
    print r.get_ruler_data()
    
    r = ruler(0, 30, (0, 100))
    print r.get_ruler_data()

    r = ruler(100, 1000, (0, 100))
    print r.get_ruler_data()

    r = ruler(10000000, 15000000, (0, 500))
    print r.get_ruler_data()
    
    r = ruler(10000000, 10200000, (0, 500))
    print r.get_ruler_data()
    