"""
gDraw, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

NOTES:
------
. This needs to be split, at the moment a lot of interface code is spilling
    over into gDraw. needs to be split into a strictly drawing backend,
    and an interface class.
. gDraw is a wrapper for whatever vector-draw backend is in use.
. The vector back-end is not normally available outside of gDraw as its
    implementation is likely to change.
. only gDraw external procedures are valid.
. procedures prefixed with _ are undocumented, internal and liable to change.
  (semi-private classes)
. DXF/PS exporter? Yes please.
. make a set of drawing primitives in case we need to swap out the back renderer
    some time in the future.
. We're currently using the 'toy' text-rendering API.
    Is this adequate for our needs?
    - will only render UTF-8 fonts... Anything else requried?
    that rules out greek symbols, Chinese, Cyrillic etc...? Right?
"""

import sys, os, math

import opt

from error import *
from constants import *
from boundbox import bBox
from operator import itemgetter
from glbase_wrapper import location

MAX_TRACKS = 10 # maximum number of tracks

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
    def __init__(self, genome):
        """
        initialise the gDraw
        pass a wx.Panel that the drawer can atach to.
        It will bind the OnPaint Event for that panel. as well.
        returns the valid panel for drawing into and binding.
        """
        self.validDraw = False # set to true if it is valid to draw
        self.scale = 1
        self.colBoxes = [] # list of boundbox objects.
        self.paintQ = []
        self.tracks = [] # currently visible tracks
        self.trackBoxes = [False for n in xrange(MAX_TRACKS)] # list of TrackBoxes in use.

        self.genome = genome

        self.chromosome = "1"
        self.lbp = 1
        self.rbp = 1000

        # set up dummy values for the view
        self.w = 100
        self.h = 200

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

    def bindTrack(self, track):
        """
        bind a drawing track extra to genome.
        """
        self.tracks.append({"data": track, "track_location": self._getNextTrackBox()})

    def _getNextTrackBox(self):
        """
        get the next available track location and return the counding coordinates of the block.
        """
        for index, track in enumerate(self.trackBoxes):
            if not track:
                self.trackBoxes[index] = True
                return(-(index * opt.track.height_px)) # 60 = genome track

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
        self.chromosome = str(chromosome)
        # test for valid chr.
        self.lbp = leftBasePair
        self.rbp = rightBasePair
        self.scale = (rightBasePair - leftBasePair) / 10000
        self.delta = self.rbp - self.lbp
        self.deltaf = float(self.delta)
        self.bps_per_pixel = self.delta / float(self.w)
        print self.bps_per_pixel, self.delta, self.w

        self.curr_loc = location(chr=self.chromosome, left=self.lbp, right=self.rbp)

        # get the new paintQ:
        self.paintQ = self.genome.getAllDrawableFeaturesInRange(self.curr_loc)
        for track in self.tracks:
            #try:
            self.paintQ.append({"type": "graph",
                "array": track["data"].get_array(location(loc=self.curr_loc), resolution=self.bps_per_pixel, read_extend=150),
                "track_location": track["track_location"],
                "name": track["data"].name})
            #except: # track probably doens't have this chr?
            #    pass
        self.paintQ.reverse()
        print [i["type"] for i in self.paintQ]
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
        loc = self._realToLocal(0, 0)
        self.ctx.move_to(loc[0], loc[1]-2)
        self.ctx.line_to(self.w, loc[1]-2)
        self.ctx.move_to(0, loc[1]+2)
        self.ctx.line_to(self.w, loc[1]+2)
        self.ctx.set_line_width(1.5)
        self.ctx.stroke()

    def _drawPoint(self, location, data):
        """
        drawFeatures of the type "Point"
        data should be a list of coordinates within the location span
        """
        pass

    def _drawSpan(self, location, data):
        """
        drawFeatures of the type "Span"
        data should be a list of coordinates of the form chrX:left-right within the location span
        """
        pass

    def drawObject(self, data):
        """
        add an object onto the Q
        """
        self.paintQ.append(data)

    def drawRuler(self):
        """
        draw the ruler.

        """
        if not self.validDraw: raise CairoDrawError

        x,y,w,h = self._getTextExtents("Chromosome %s" % str(self.chromosome))
        self._drawText(5, opt.ruler.height_px + 22, opt.ruler.font, "Chromosome %s" % str(self.chromosome), size=12)

        self._setPenColour(opt.ruler.colour)
        # work out a good scale representation
        # wether to draw at 100, 1000, 10000, 100000, 1000000 ...

        # current view delta = self.delta

        a = round(self.delta, 1) # get the nearest 1XXXXX .. XXX

        # ten thousands
        for index, window_size in enumerate([a/100, a/10, a]):

            nearest = int(math.ceil(float(self.lbp+1) / window_size) * window_size)
            self.ctx.set_line_width(opt.ruler.line_width * index+0.5)

            for real_offset in xrange(nearest, self.rbp, int(window_size)):
                screen_offset = (self.w * (float(real_offset - self.lbp) / self.delta))
                self.ctx.move_to(screen_offset, 0)
                self.ctx.line_to(screen_offset, opt.ruler.height_px * index+0.5)
                self.ctx.stroke()
                if index == 1: # write numbers at 1/10 scale.
                    self._drawText(screen_offset +2, opt.ruler.text_height, opt.ruler.font, str(real_offset), opt.ruler.font_size)

        # set up the bounding box.
        self.colBoxes.append(bBox((0,0,self.w, opt.ruler.text_height), (0,0), None, None))
        return(True)

    def exportPostScript(self, _filename):
        pass

    def exportDXF(self, _filename):
        pass

    def exportPNG(self, _filename):
        pass

    def exportBMP(self, _filename):
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
    #------------------------------------------------------------------

    def _realToLocal(self, x, y):
        """
        (Internal)
        Convert real genomic coords to local pixel coordinates.
        local to the location of the genome line.
        """
        return((self.w * ((x-self.lbp) / self.deltaf), (self.h - 30) + y))

    def _localToReal(self, sx, sy):
        """
        opposite of realToLocal()
        """
        x = 0
        y = 0
        return((x,y))

    def _getTextExtents(self, text):
        """
        (Internal)
        returns a bounding box of the text return = (x,y,w,h)
        """
        return(self.ctx.text_extents(text)[:4])

    def _drawTrackBackground(self, track_location):
        """
        track_location is the bottom edge of the track block
        """
        # get an available track slot
        self._setPenColour( (0.95,0.95,0.95) )
        base_loc = self._realToLocal(0, track_location)
        self.ctx.rectangle(0, base_loc[1]-opt.track.height_px, self.w, opt.track.height_px-2) # 30 = half genomic track size
        self.ctx.fill()

    def _drawGraphTrack(self, track_data):
        # need to resample the array to speed up the draw
        self._drawTrackBackground(track_data["track_location"])
        self._setPenColour( (0,0,0) )
        self.ctx.set_line_width(0.5)
        coords = []
        lastpx = -1
        for index, value in enumerate(track_data["array"]):
            loc = self._realToLocal(self.lbp + index, track_data["track_location"])
            #if int(loc[0]) > lastpx: # this means only draw one per pixel
            #    lastpx = loc[0]
            #    coords.append( (index, loc[1] - 30 - value)) # +30 locks it to the base of the track
            coords.append( (index, loc[1] - value)) # +30 locks it to the base of the track

        self.ctx.move_to(coords[0][0], coords[0][1]) # start x,y
        for index, item in enumerate(coords):
            self.ctx.line_to(item[0], item[1])
        self.ctx.stroke()
        self._drawText(0, loc[1] - 15 , opt.graphics.font, track_data["name"])

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

    def _drawGene(self, data):
        """
        draw Features of the type "Gene"
        should be an dict containing the following keys:
        type: gene
        loc: location span of gene
        strand: strand of gene
        cds_loc: cds loc span
        exonStarts: list of exon start locations
        exonEnds: list of ends
        """
        if not self.validDraw: raise ErrorCairoDraw

        #print "t:", ((data["left"]-self.lbp) / self.deltaf), ((data["right"]-self.lbp) / self.deltaf)

        posLeft = self._realToLocal(data["loc"]["left"], 0)[0]
        posRight = self._realToLocal(data["loc"]["right"], 0)[0]
        self.ctx.set_line_width(0.5)
        self._setPenColour(opt.graphics.gene_colour)

        #---------------------------------------------------------------
        # draw gene blocks.
        tc = [] # build a list of the genome coordinates.
        for item in data["exonStarts"]:
            tc.append({"c": item ,"t": "es"})
        for item in data["exonEnds"]:
            tc.append({"c": item, "t": "ee"})
        tc.append({"c": data["cds_loc"]["left"], "t": "cdss"})
        tc.append({"c": data["cds_loc"]["right"], "t": "cdse"})

        tc = sorted(tc, key=itemgetter("c"))

        current_offset = opt.graphics.gene_height
        coords = []
        coords.append(self._realToLocal(data["loc"]["left"], 0))
        coords.append(self._realToLocal(data["loc"]["left"], + opt.graphics.gene_height))
        for c in tc:
            if c["t"] == "cdss":
                current_offset = opt.graphics.cds_height
                coords.append(self._realToLocal(c["c"], + opt.graphics.gene_height))
                coords.append(self._realToLocal(c["c"], + opt.graphics.cds_height))
            elif c["t"] == "cdse":
                current_offset = opt.graphics.gene_height
                coords.append(self._realToLocal(c["c"], + opt.graphics.cds_height))
                coords.append(self._realToLocal(c["c"], + opt.graphics.gene_height))
            elif c["t"] == "es":
                coords.append(self._realToLocal(c["c"], 0))
                coords.append(self._realToLocal(c["c"], + current_offset))
            elif c["t"] == "ee":
                coords.append(self._realToLocal(c["c"], + current_offset))
                coords.append(self._realToLocal(c["c"], 0))

        coords.append(self._realToLocal(data["loc"]["right"], + opt.graphics.gene_height))
        coords.append(self._realToLocal(data["loc"]["right"], 0))

        self.ctx.move_to(coords[0][0], coords[0][1])
        for index, item in enumerate(coords):
            self.ctx.line_to(item[0], item[1])
        self.ctx.move_to(coords[0][0], coords[0][1])

        #coords.reverse()
        for item in coords:
            self.ctx.line_to(item[0], self.h - 30 - (item[1] - (self.h-30)))
        #self.ctx.stroke()
        self.ctx.fill()

        self.ctx.set_line_width(1)

        #---------------------------------------------------------------
        # Draw gene arrow

        if data["strand"] == "+": # top strand
            loc = self._realToLocal(data["loc"]["left"], 0)
            # arrow.
            self.ctx.move_to(loc[0]+10, loc[1]-20)
            self.ctx.line_to(loc[0]+10, loc[1]-30)
            self.ctx.line_to(loc[0]+20, loc[1]-20)
            self.ctx.line_to(loc[0]+10, loc[1]-10)
            self.ctx.line_to(loc[0]+10, loc[1]-20)
            self.ctx.stroke()
            self._drawText(loc[0], loc[1]-25, "Arial", data["name"], 12)
        elif data["strand"] == "-":
            loc = self._realToLocal(data["loc"]["right"], 0)
            # arrow.
            self.ctx.move_to(loc[0]-10, loc[1]+20)
            self.ctx.line_to(loc[0]-10, loc[1]+30)
            self.ctx.line_to(loc[0]-20, loc[1]+20)
            self.ctx.line_to(loc[0]-10, loc[1]+10)
            self.ctx.line_to(loc[0]-10, loc[1]+20)
            self.ctx.stroke()
            self._drawText(loc[0], loc[1]+25, "Arial", data["name"], 12)
        else:
            raise ErrorInvalidGeneDefinition

        return(True)

    def _paint(self, ctx):
        """
        (Internal)
        painter procedure. pass a valid Cairo context for drawing onto.

        draws every object in paintQ
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

        func_dict = {
            "gene": self._drawGene,
            "lncRNA": self._drawGene,
            "microRNA": self._drawGene,
            "graph": self._drawGraphTrack
        }

        self.drawChr(None)

        for item in self.paintQ:
            if func_dict.has_key(item["type"]):
                func_dict[item["type"]](item)

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

    def isColliding(self, x, y):
        """
        returns a collision type (see boundbox)
        or FALSE
        """
        for box in self.colBoxes:
            c = box.collideB()
            if c:
                print c["type"]
                return(c["type"])
