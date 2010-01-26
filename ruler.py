"""

ruler.py

abstraction of a generic ruler class.

"""

from __future__ import division

class ruler:
    def __init__(self, min, max, display_pixel_size, useSI_label=True):
        """
        **Purpose**
            an abstract class to represent rulers.

            (minimum values not implemented yet)

        **Arguments**
            min

            max

            display_pixel_size

            useSI_label (default = True)
                append k for 1000, M for 1,000,000, G for 1,000,000,000 etc...
                to the labels.

        **Results**
            A valid ruler object
        """
        self.min = min
        self.max = max
        self.display_px = display_pixel_size
        self.__recalculate_fractions()

    def __recalculate_fractions(self):
        """
        (Internal)
        recalculate things like units per pixel and other
        internally used values.
        """
        #self.__px_in_units = self.max * (self. # the number of pixels per unit
        self.delta = self.max - self.min
        self.__units_in_px = self.max / self.display_px # number of units per pixel

    def set(self, min, max):
        """
        **Purpse**

        **Arguments**
            min
                update the minimum value

            max
                update the maximum value
        """
        self.min = min
        self.__max = maximum_value
        self.__recalculate_fractions()

    def get_ruler_data(self, percent=20):
        """
        **Purpose**

        **Arguments**

            percent (Optional)
                the percent to use to break up the ruler.
                (defaults to 20%)

        **Returns**
        sends back a list containing tuple pairs in the form:
        [(location_px, "label"), (..) .. (..)]
        """
        fp = percent / 100.0

        ppx = round(self.display_px * fp) # work out how many pixels 20% is.
        p = round(self.max * fp) # work out how many pixels 20% is.

        print p, ppx

        best = 0
        minim = 100000000000

        # find the closest unit to p:
        scales_and_labels = {
            0.001: "m",
            0.01: "m",
            0.1: "m",
            1: "",
            10: "",
            100: "",
            1000: "k",
            10000: "k",
            100000: "k",
            1000000: "M",
            10000000: "M",
            100000000: "M",
            1000000000: "G",
            10000000000: "G",
            100000000000: "G", # don't support more than this.
            1000000000000: "?" # what's bigger than gig?
            }
            # a 64-bit maxint:
            #9223372036854775807

        for i, v in enumerate(scales_and_labels):
            if minim > abs(p - v):
                minim = abs(p - v)
                best = v

        #print best, scales_and_labels[best]

        res = []
        for v in xrange(0, self.max, best):# min vals not supported, floats will break this?
            px_location = self.min + (self.display_px * (v / self.delta))
            res.append( (px_location, "%s%s" % (v, scales_and_labels[best])) )

        return(res)

    def draw(self, cairo_context, draw_top_left_tuple):
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
        draw_stuff = self.get_ruler_data()

        x,y,w,h = self.__getTextExtents("Chromosome %s" % str(self.chromosome))
        self.__drawText(5, opt.ruler.height_px + 22, opt.ruler.font, "Chromosome %s" % str(self.chromosome), size=12)

        self.__setPenColour(opt.ruler.colour)
        # work out a good scale representation
        # wether to draw at 100, 1000, 10000, 100000, 1000000 ...

        # current view delta = self.delta

        a = round(self.delta, 1) # get the nearest 1XXXXX .. XXX

        # ten thousands
        for index, window_size in enumerate([int(a/100), int(a/10), int(a)]):

            nearest = int(math.ceil(float(self.lbp+1) / window_size) * window_size)
            self.ctx.set_line_width(opt.ruler.line_width * index+0.5)

            for real_offset in xrange(nearest, self.rbp, int(window_size)):
                screen_offset = (self.w * (float(real_offset - self.lbp) / self.delta))
                self.ctx.move_to(screen_offset, 0)
                self.ctx.line_to(screen_offset, opt.ruler.height_px * index+0.5)
                self.ctx.stroke()
                #if index == 1: # write numbers at 1/10 scale.
                #    self.__drawText(screen_offset +2, opt.ruler.text_height, opt.ruler.font, str(real_offset), opt.ruler.font_size)

        return((0,0,self.w, opt.ruler.text_height)) # return the colbox
        return(None)

if __name__ == "__main__":
    r = ruler(5, 10, 100)
    print r.get_ruler_data()

    r = ruler(0, 30, 100)
    print r.get_ruler_data()

    r = ruler(100, 1000, 100)
    print r.get_ruler_data()
