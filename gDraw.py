"""
gDraw, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

NOTES:
------
. gDraw is a wrapper for whatever vector-draw backend is in use.
. The vector back-end is not normally available outside of gDraw as its
    implementation is likely to change.
. only gDraw external procedures are valid.
. procedures prefixed with _ are undocumented, internal and liable to change.
  (private classes)
. DXF/PS exporter? Yes please.
. make a set of drawing primitives in case we need to swap out the back renderer
    some time in the future.
. We're currently using the 'toy' text-rendering API.
    Is this adequate for our needs?
    - will only render UTF-8 fonts... Anything else requried?
    that rules out greek symbols, Chinese Cyrillic etc...?
"""

import sys, os, math, opt

from error import *
from constants import *
from boundbox import bBox

#----------------------------------------------------------------------
# The Cairo interface.
import wx.lib.wxcairo

import ctypes
import cairo
from ctypes.util import find_library

cairo_dll = ctypes.CDLL(find_library("cairo"))

# Pycairo's API representation (from pycairo.h)
class Pycairo_CAPI(ctypes.Structure):
   _fields_ = [
      ('Context_Type', ctypes.py_object),
      ('Context_FromContext', ctypes.PYFUNCTYPE(ctypes.py_object,
                                                ctypes.c_void_p,
                                                ctypes.py_object,
                                                ctypes.py_object)),
      ('FontFace_Type', ctypes.py_object),
      ('FontFace_FromFontFace', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('FontOptions_Type', ctypes.py_object),
      ('FontOptions_FromFontOptions', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Matrix_Type', ctypes.py_object),
      ('Matrix_FromMatrix', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Path_Type', ctypes.py_object),
      ('Path_FromPath', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Pattern_Type', ctypes.py_object),
      ('SolidPattern_Type', ctypes.py_object),
      ('SurfacePattern_Type', ctypes.py_object),
      ('Gradient_Type', ctypes.py_object),
      ('LinearGradient_Type', ctypes.py_object),
      ('RadialGradient_Type', ctypes.py_object),
      ('Pattern_FromPattern', ctypes.c_void_p),
      ('ScaledFont_Type', ctypes.py_object),
      ('ScaledFont_FromScaledFont', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Surface_Type', ctypes.py_object),
      ('ImageSurface_Type', ctypes.py_object),
      ('PDFSurface_Type', ctypes.py_object),
      ('PSSurface_Type', ctypes.py_object),
      ('SVGSurface_Type', ctypes.py_object),
      ('Win32Surface_Type', ctypes.py_object),
      ('XlibSurface_Type', ctypes.py_object),
      ('Surface_FromSurface', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Check_Status', ctypes.PYFUNCTYPE(ctypes.c_int, ctypes.c_int))]

# look up the API
ctypes.pythonapi.PyCObject_Import.restype = ctypes.POINTER(Pycairo_CAPI)
pycairo_api = ctypes.pythonapi.PyCObject_Import("cairo", "CAPI").contents;

ContextType = pycairo_api.Context_Type

def Context_FromSWIGObject(swigObj):
    """
    (Internal)
    get the Cairo object for drawing.
    """
    ptr = ctypes.c_void_p(int(swigObj))
    #increment the native context's ref count, since the Pycairo_Context decrements it
    #when it is finalised.
    cairo_dll.cairo_reference(ptr)
    return pycairo_api.Context_FromContext(ptr, ContextType, None)

class drawPanel(wx.Panel):
    """
    An empty Panel class that sets things up suitbly for Cairo.
    """
    def __init__(self, parent, paint_Procedure):
        wx.Panel.__init__(self, parent, -1, size=(-1,-1), style=wx.FULL_REPAINT_ON_RESIZE) # This has to be set to full repaint :)
        self.Bind(wx.EVT_PAINT, paint_Procedure, self)
        self.Fit()

#----------------------------------------------------------------------
# End of the Cairo interface

class gDraw:
    def __init__(self):
        """
        initialise the gDraw
        pass a wx.Panel that the drawer can atach to.
        It will bind the OnPaint Event for that panel. as well.
        returns the valid panel for drawing into and binding.
        """
        self.validDraw = False # set to true if it is valid to draw
        self.scale = 1
        self.colBoxes = [] # list of boundbox objects.

    def _debug_draw_col_boxes(self):
        """
        (Internal)
        Debug routine to draw the locations of the collision boxes.
        """
        self._setPenColour((0.6, 0.2, 0, 0.5))
        for box in self.colBoxes:
            dim = box.getDimensions()
            self.ctx.rectangle(dim[0], dim[1], dim[2], dim[3])
            self.ctx.fill()
        self.ctx.stroke()

    def bindPanel(self, panel):
        """
        bind a wx.Panel into gDraw.
        """
        self.panel = drawPanel(panel, self.OnPaint)
        return(self.panel)

    def setViewPortSize(self, w, h):
        """
        set the size of the viewport;
        """
        self.w = w
        self.h = h
        self.aspect = abs(float(w) / h)
        self.halfw = w / 2
        self.halfh = h / 2

    def setLocation(self, chromosome, leftBasePair, rightBasePair):
        """
        Set the location of the view, this will show the left most and rightmost
        base pair number.
        It also calculates the internal scale representations.
        """
        try:
            self.chromosome = int(chromosome)
        except ValueError:
            if chromosome in validChrNames: # from constants
                self.chromosome = chromosome
            else:
                raise ErrorInvalidChromosome

        self.lbp = leftBasePair
        self.rbp = rightBasePair
        self.scale = (rightBasePair - leftBasePair) / 10000
        self.delta = self.rbp - self.lbp
        self.deltaf = float(self.delta)
        return(True)

    def setDrawAttribute(self, attribute, value):
        """
        setDrawAttribute
        """
        pass

    def drawChr(self, location_span):
        """
        draw the basic chromosome
        """
        if not self.validDraw: raise CairoDrawError
        self.ctx.set_source_rgb(0, 0, 0)
        self.ctx.move_to(0, (self.h / 2)-1)
        self.ctx.line_to(self.w, (self.h / 2)-1)
        self.ctx.move_to(0, (self.h / 2)+2)
        self.ctx.line_to(self.w, (self.h / 2)+2)
        self.ctx.set_line_width(2)
        self.ctx.stroke()

    def drawPoint(self, location, data):
        """
        drawFeatures of the type "Point"
        data should be a list of coordinates within the location span
        """
        pass

    def drawSpan(self, location, data):
        """
        drawFeatures of the type "Span"
        data should be a list of coordinates of the form chrX:left-right within the location span
        """
        pass

    def drawRuler(self):
        """
        draw the ruler.

        This needs to inteligently scale based on the current number of bp's displayed.
        """
        if not self.validDraw: raise CairoDrawError

        x,y,w,h = self._getTextExtents("Chromosome %s" % str(self.chromosome))
        self._drawText(5, opt.ruler.text_height+w, opt.ruler.font, "Chromosome %s" % str(self.chromosome))

        self._setPenColour(opt.ruler.colour)
        # work out a good scale representation
        # wether to draw at 100, 1000, 10000, 100000, 1000000 ...

        # ten housands
        for index, window_size in enumerate([100, 1000, 10000]):

            nearest = int(math.ceil(float(self.lbp+1) / window_size) * window_size)
            self.ctx.set_line_width(opt.ruler.line_width * index+0.5)

            for real_offset in xrange(nearest, self.rbp, window_size):
                screen_offset = (self.w * (float(real_offset - self.lbp) / self.delta))
                self.ctx.move_to(screen_offset, 0)
                self.ctx.line_to(screen_offset, opt.ruler.height_px * index+0.5)
                self.ctx.stroke()
                if index == 1:
                    self._drawText(screen_offset +2, opt.ruler.text_height, opt.ruler.font, str(real_offset), opt.ruler.font_size)

        # set up the bounding box.
        self.colBoxes.append(bBox((0,0,self.w, opt.ruler.text_height), (0,0), None, None))
        return(True)

    def drawGene(self, data):
        """
        drawFeatures of the type "Gene"
        data should be a dict of the form:
        [chr, left, right,
        exon_coords(local), strand, cds_start, cds_end]
        """
        if not self.validDraw: raise ErrorCairoDraw

        # if the gene is not within the viewport the quit out

        #print "t:", ((data["left"]-self.lbp) / self.deltaf), ((data["right"]-self.lbp) / self.deltaf)

        # should turn this pos stuff into a macro
        posLeft = self.w * ((data["left"]-self.lbp) / self.deltaf)
        posRight = (self.w * ((data["right"]-self.lbp) / self.deltaf))
        #print posLeft, posRight, self.delta

        self.ctx.set_line_width(1)

        if data["type"] == "gene":
            self._setPenColour(opt.graphics.gene_colour)
            hei_off = opt.graphics.gene_height # height offset
        elif data["type"] == "lncRNA":
            self._setPenColour(opt.graphics.lncrna_colour)
            hei_off = opt.graphics.lncrna_height
        elif data["type"] == "microRNA":
            self._setPenColour(opt.graphics.microRNA_colour)
            hei_off = opt.graphics.microRNA_height
        else:
            self.ctx.set_source_rgb(0,0,0)
            hei_off = 5

        if data["strand"] == "+": # top strand
            self.ctx.rectangle(posLeft,self.halfh-hei_off,posRight - posLeft,hei_off * 2)
            self.ctx.fill()
            # the line, up.
            self.ctx.move_to(posLeft, self.halfh)
            self.ctx.line_to(posLeft, self.halfh-20)
            self.ctx.line_to(posLeft+20, self.halfh-20)
            self.ctx.stroke()
            # arrow.
            self.ctx.move_to(posLeft+10, self.halfh-20)
            self.ctx.line_to(posLeft+10, self.halfh-30)
            self.ctx.line_to(posLeft+20, self.halfh-20)
            self.ctx.line_to(posLeft+10, self.halfh-10)
            self.ctx.line_to(posLeft+10, self.halfh-20)
            self.ctx.stroke()
        elif data["strand"] == "-":
            self.ctx.rectangle(posLeft,self.halfh-hei_off,posRight - posLeft,hei_off * 2)
            self.ctx.fill()
            # the line
            self.ctx.move_to(posRight, self.halfh)
            self.ctx.line_to(posRight, self.halfh+20)
            self.ctx.line_to(posRight-20, self.halfh+20)
            self.ctx.stroke()
            # arrow.
            self.ctx.move_to(posRight-10, self.halfh+20)
            self.ctx.line_to(posRight-10, self.halfh+30)
            self.ctx.line_to(posRight-20, self.halfh+20)
            self.ctx.line_to(posRight-10, self.halfh+10)
            self.ctx.line_to(posRight-10, self.halfh+20)
            self.ctx.stroke()
        else:
            raise ErrorInvalidGeneDefinition
            # bottom strand

    def exportPostScript(self, _filename):
        pass

    def exportDXF(self, _filename):
        pass

    def exportPNG(self, _filename):
        pass

    def exportBMP(self, _filename):
        pass

    def exportJPG(self, _filename):
        pass

    def getPanel(self):
        """
        get the current panel.
        """
        return(self.panel)

    def forceRedraw(self):
        """
        Use me to force a repaint, rather than calling OnPaint directly.
        """
        self.OnPaint(None)
        return(True)

    def OnPaint(self, event):
        """
        Call to get Cairo to paint, this is done automatically if bound
        to a Panel.
        """
        # override;
        try:
            dc = wx.PaintDC(self.panel) # make each time OnPaint is called.
            #dc = wx.AutoBufferedPaintDC(self) # not clear why this doesn't work ...
            self.size = dc.GetSizeTuple()
            self.setViewPortSize(self.size[0], self.size[1])
            gc = wx.GraphicsContext.Create(dc)
            nc = gc.GetNativeContext()
            ctx = Context_FromSWIGObject(nc)
        except:
            raise ErrorCairoAcquireDevice

        self._paint(ctx)
        return(True)

    #------------------------------------------------------------------
    # Internal painters
    def _getTextExtents(self, text):
        """
        (Internal)
        returns a bounding box of the text return = (x,y,w,h)
        """
        return(self.ctx.text_extents(text)[:4])

    def _drawText(self, x, y, font, text, size=12, colour=(0,0,0), style=None):
        """
        (Internal - helper)
        Draw text to the screen, font is a string for the font name.
        Does not check for availability of the font.
        Will accept several styles, see constants for details.

        """
        self._setPenColour(colour)
        self.ctx.select_font_face(font, txtToCairoA[style], txtToCairoB[style])
        self.ctx.move_to(x, y)
        self.ctx.show_text(text)

    def _paint(self, ctx):
        """
        (Internal)
        painter procedure. pass a valid Cairo context for drawing onto.
        """
        if ctx:
            self.validDraw = True
        else:
            raise ErrorCairoAcquireDevice
        self.ctx = ctx

        # kill all the colission boxes
        self.colBoxes = []

        # blank the screen:
        self._setPenColour(opt.graphics.screen_colour)
        ctx.rectangle(0,0,self.w,self.h)
        ctx.fill()

        # do the drawing's here:
        #self.setLocation(3, 9500, 15000) # what the viewport can see.
        #self.setLocation(3, 5000, 30000)
        self.drawChr(None)
        tgene = {"type": "gene", "strand": "+", "left": 12344, "right": 13455}
        self.drawGene(tgene)
        tgene = {"type": "lncRNA", "strand": "-", "left": 10000, "right": 11000}
        self.drawGene(tgene)
        tgene = {"type": "microRNA", "strand": "+", "left": 11050, "right": 11100}
        self.drawGene(tgene)

        self.drawRuler()

        # render and finish drawing.
        self.validDraw = False

        # any further (internal) drawing goes here.
        if opt.debug.draw_collision_boxes: self._debug_draw_col_boxes()
        return(True)

    def _setPenColour(self, colour):
        """
        (Internal - draw primitive)
        A macro for changing the current pen colour.
        deals with rgb and rgba intelligently.
        """
        if len(colour) == 3:
            self.ctx.set_source_rgb(colour[0], colour[1], colour[2])
        elif len(colour) == 4:
            self.ctx.set_source_rgba(colour[0], colour[1], colour[2], colour[3])
        else:
            return(False)
        return(True)
