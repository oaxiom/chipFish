"""
gDraw, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

NOTES:
------
. This needs to be split into two
    at the moment a lot of interface code is spilling
    over into gDraw. needs to be split into a strictly drawing backend,
    and an interface class. Actually, why? Cairo is mature...

. I think this needs to be split further than 2.
    gDraw should conatain the gui draw interface, and the actual drawing primitives __draw* moved into a
    _draw module.

. gDraw is the only location to use cairo... Actually, why restrict this? gDraw is getting far too big...
. Currently using the 'toy' text-rendering API.
    Is this adequate for our needs?
    - will only render UTF-8 fonts... Anything else requried?
    that rules out greek symbols, Chinese, Cyrillic etc...? Right?

TODO:
-----
. change the location to a <location> and use that datatype instead (cleaner)
. DXF/PS exporter?
"""

from __future__ import division

import sys, os, math

import opt

from error import *
from constants import *
from boundbox import bbox
from operator import itemgetter
from glbase_wrapper import location
from data import *
from ruler import ruler

MAX_TRACKS = 10 # maximum number of tracks - should be in config

import wx.lib.wxcairo
import cairo

class drawPanel(wx.Panel):
    """
    An empty Panel class that sets things up suitably for Cairo.
    """
    def __init__(self, parent, paint_Procedure):
        wx.Panel.__init__(self, parent, -1, size=(-1,-1), style=wx.FULL_REPAINT_ON_RESIZE) # This has to be set to full repaint :)
        self.Bind(wx.EVT_PAINT, paint_Procedure, self)
        self.Fit()

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

        # set up dummy values for the view these get updated to sensible values later
        self.w = 100
        self.h = 200

    def OnPaint(self, event, cairo_context=None):
        """
        **Event**
            Call to get Cairo to paint, this is done automatically if bound
            to a Panel. You do not need to call this explicitly.
            If you pass a Cairo Context cairo will paint to that
            rather than generate a context for the gui.
        """
        # try to get cairo:
        if not cairo_context:
            try:
                dc = wx.PaintDC(self.panel) # make each time OnPaint is called.
                #dc = wx.AutoBufferedPaintDC(self) # not clear why this doesn't work ...
                self.size = dc.GetSizeTuple()
                self.setViewPortSize(self.size[0], self.size[1])
                cairo_context = wx.lib.wxcairo.ContextFromDC(dc)
            except:
                raise ErrorCairoAcquireDevice

        self.ctx = cairo_context

        draw_modes_dict = {
            "gene": self.__drawGene,
            "lncRNA": self.__drawGene,
            "microRNA": self.__drawGene,
            "graph": self.__drawTrackGraph,
            "graph_split_strand": self.__drawTrackGraph_split_strand,
            "bar": self.__drawTrackBar,
            "spot": self.__drawTrackSpot
            }

        self.delta = self.rbp - self.lbp
        self.deltaf = float(self.delta)
        self.bps_per_pixel = self.delta / float(self.w)

        self.curr_loc = location(chr=self.chromosome, left=self.lbp, right=self.rbp)

        # kill all the colission boxes
        self.__col_boxs = []

        # blank the screen:
        self.__setPenColour(opt.graphics.screen_colour)
        self.ctx.rectangle(0,0,self.fullw,self.h)
        self.ctx.fill()

        if opt.draw.double_lines_for_genome:
            self.__drawChr(None)

        # get the new paintQ:
        # draw the genome:
        genome_items = self.genome.getAllDrawableFeaturesInRange(self.curr_loc)
        for item in genome_items:
            draw_modes_dict[item["type"]](item)

        # collect the data for the tracks and draw them on the screen.
        for track in self.tracks:
            # basic data:
            draw_data = {"type": track["type"],
                "track_location": track["track_location"],
                "name": track["data"].name}

            if track["type"] == "graph":
                draw_data["array"] = track["data"].get_data("graph", location(loc=self.curr_loc),
                    resolution=self.bps_per_pixel, read_extend=50)
            elif track["type"] == "bar":
                draw_data["array"] = track["data"].get_data("bar", location(loc=self.curr_loc),
                    resolution=self.bps_per_pixel, read_extend=0)
            elif track["type"] == "spot":
                draw_data["array"] = track["data"].get_data("spot", location(loc=self.curr_loc))
            elif track["type"] == "graph_split_strand":
                draw_data["array"] = track["data"].get_data("graph", location(loc=self.curr_loc),
                    strand=True, resolution=self.bps_per_pixel, read_extend=0)

            # and draw:
            colbox = draw_modes_dict[draw_data["type"]](draw_data)

            # add a collision item if
            if colbox:
                self.__col_boxs.append(bbox(colbox, track, "track"))

        self.__col_boxs.append(bbox(self.__drawRuler(), None, "ruler"))
        if opt.draw.scale_bar:
            self.__drawScaleBar()

        # any further (internal) drawing goes here.
        if opt.debug.draw_collision_boxes: self.__debug_draw_col_boxes()
        return(True)

    def __debug_draw_col_boxes(self):
        """
        (Internal)
        Debug routine to draw the locations of the collision boxes.

        set me in opt.debug
        """
        self.__setPenColour((0.6, 0.2, 0, 0.3))
        for box in self.__col_boxs:
            self.ctx.rectangle(*box.get_dimensions())
            self.ctx.fill()
        self.ctx.stroke()

    def bindPanel(self, panel):
        """
        This is special usage, and not so intuitive at the moment...
        Maybe in the future it will be better.
        bind a wx.Panel into gDraw.
        """
        self.panel = drawPanel(panel, self.OnPaint)
        return(self.panel)

    def bindTrack(self, track, track_type=None):
        """
        bind a drawing track extra to genome.

        # valid track types:
        graph
        bar
        spot
        """
        # if no track_type try to guess from the track object
        if not track_type:
            track_type = track._default_draw_type

        if track_type not in valid_track_draw_types:
            raise ErrorTrackDrawTypeNotFound, track_type

        self.tracks.append({"data": track, "track_location": self.__getNextTrackBox(track_type), "type": track_type})

    def __getNextTrackBox(self, track_type):
        """
        get the next available track location and return the counding coordinates of the block.
        """
        currentLoc = 0
        for index, track in enumerate(self.trackBoxes):
            if track:
                currentLoc += opt.track.height_px[track]
            if not track:
                self.trackBoxes[index] = track_type
                return(-(index + (currentLoc))-opt.track.genome_base_offset) # 60 = genome track

    def setViewPortSize(self, w, h):
        """
        set the size of the viewport;
        """
        self.fullw = w # for blanking screen
        self.w = w - opt.graphics.right_border_width # a small right most border for editing trakcs
        self.h = h
        self.aspect = abs(float(self.w) / self.h)
        self.halfw = w / 2
        self.halfh = h / 2

    def setLocation(self, chromosome=None, leftBasePair=None, rightBasePair=None, loc=None, **kargs):
        """
        **Purpose**
            Set the location of the view, this will show the left most and rightmost
            base pair number.
            It also calculates the internal scale representations.

        **Arguments**
            This method is a little scizophrenic at the moment, supporting both
            an old-style and new-style <location> based associaton
            later it should only support the new-style <location>

        **Returns**
            Nothing
            Does not rebuild the Cairo Display! but it makes the next call to
            OnPaint() or forceRedraw() correct.
        """
        if chromosome: # old-style assignation
            self.chromosome = str(chromosome)
            self.lbp = leftBasePair
            self.rbp = rightBasePair
        elif loc: # new-style <location> assignation.
            self.chromosome = loc["chr"]
            self.lbp = loc["left"]
            self.rbp = loc["right"]
        # sanity checking? Neccesary?

    def getLocation(self):
        """
        **Purpose**
            get the current view genome location

        **Arguments**
            None

        **Returns**
            returns a <location>
        """
        return(location(chr=self.chromosome, left=self.lbp, right=self.rbp))

    def move(self, mode="centre", percent=5):
        """
        move the view by a fixed percent according to 'mode'
        valid modes are:
            left = left
            right = right
            zoomout = zoom out
            zoomin = zoom in
        will move the display by n percent.
        """
        if mode not in valid_move_modes:
            return(False)

        # get the current bp move percent.
        move_percent = int((self.rbp - self.lbp) * (percent / 100.0))

        if mode == "right":
            self.lbp += move_percent
            self.rbp += move_percent
        elif mode == "left":
            self.lbp -= move_percent
            self.rbp -= move_percent
        elif mode == "zoomout":
            self.lbp -= move_percent
            self.rbp += move_percent
        elif mode == "zoomin":
            self.lbp += move_percent
            self.rbp -= move_percent
            if self.lbp >= self.rbp: # check not zoomed in too far.
                mid_point = self.rbp + self.lbp / 2
                self.lbp = mid_point - 1
                self.rbp = mid_point + 1
        return(True)

    def __drawChr(self, location_span):
        """
        draw the basic chromosome
        """
        self.ctx.set_source_rgb(0, 0, 0)
        loc = self.__realToLocal(0, 0)
        self.ctx.move_to(loc[0], loc[1]-2)
        self.ctx.line_to(self.w, loc[1]-2)
        self.ctx.move_to(0, loc[1]+2)
        self.ctx.line_to(self.w, loc[1]+2)
        self.ctx.set_line_width(1.5)
        self.ctx.stroke()

    def __drawScaleBar(self):
        """
        draw a scale bar to the next 000'th
        """
        # work out most relevant n'000
        # how many bp does 20% of the display take up?
        percent = round(self.delta * 0.2)
        # find the closest 10, 100, 1000 ...

        # guess the best scale describing 20% of the view
        best = 0
        res = []
        minim = 1000000000

        scales_and_labels = {10: "10bp",
            100: "100bp",
            1000: "1kbp",
            10000: "10kbp",
            100000: "100kbp",
            1000000: "1Mbp",
            10000000: "10Mbp",
            100000000: "100Mbp", # Human chr 1 is 250Mbp.
            1000000000: "1000Mbp" # Just in case there are some wierd genomes.
            }

        for i, v in enumerate(scales_and_labels):
            if minim > abs(percent - v):
                minim = abs(percent - v)
                best = v

        posLeft = self.__realToLocal(self.rbp-best, -40)
        posRight = self.__realToLocal(self.rbp, -40)

        self.ctx.set_line_width(0.5)
        self.ctx.move_to(posLeft[0]-20, posLeft[1]) # move 20px arbitrarily left
        self.ctx.line_to(posRight[0]-20, posRight[1])
        self.ctx.stroke()

        self.ctx.move_to(posLeft[0]-20, posLeft[1]-4) # move 20px arbitrarily left
        self.ctx.line_to(posLeft[0]-20, posLeft[1]+4)
        self.ctx.stroke()

        self.ctx.move_to(posRight[0]-20, posRight[1]-4) # move 20px arbitrarily left
        self.ctx.line_to(posRight[0]-20, posRight[1]+4)
        self.ctx.stroke()

        self.__drawText(posLeft[0], posRight[1]-5, opt.graphics.font, scales_and_labels[best])

    def __drawRuler(self):
        """
        draw the ruler.

        """
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
                if index == 1: # write numbers at 1/10 scale.
                    self.__drawText(screen_offset +2, opt.ruler.text_height, opt.ruler.font, str(real_offset), opt.ruler.font_size)

        return((0,0,self.w, opt.ruler.text_height)) # return the colbox

    def exportImage(self, filename, type=None):
        """
        **Purpose**

            export the current Cairo image as a 'type'
            the image type will be guessed from the filename.
            If no obvious extension is given then the string value in 'types'
            will be used.
            Finally, if that doens't make sense then it will default to a png

        **Arguments**

            filename (Required)
                filename and path (and extension) to the file to save

            type (Optional)
                the type of file to save as, this will override the
                extension in the filename. If not type is given and
                the extension doesn't make any sense then a png will be
                used.

        **Result**
            returns the actual filename used to save.
            and a file saved in filename
        """
        # guess the maximum height required:
        guess_height = abs(self.tracks[-1]["track_location"]) + opt.track.height_px[self.tracks[-1]["type"]] + opt.ruler.height_px + 30 + 32 # 32 is the 'chromosome %s' padding, 30 is some other padding I'm not really certain where it comes from...

        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 1024, guess_height)
        ctx = cairo.Context(surface)

        # preserve view variables
        oldw = self.w
        oldh = self.h
        oldfullw = self.fullw

        self.w = 1024
        self.fullw = 1024 # no bar side buttons
        self.h = guess_height

        # forceRedraw onto my surface.
        self.OnPaint(None, ctx)

        # save image
        actual_filename = filename
        filehandle = open(actual_filename, "wb")
        surface.write_to_png(filehandle)
        filehandle.close()

        # put back the gDraw view values.
        self.w = oldw
        self.h = oldh
        self.fullw = oldfullw

        # clean ups
        del ctx
        del surface
        return(actual_filename)

    def getPanel(self):
        """
        get the current panel.
        """
        return(self.panel)

    def forceRedraw(self):
        """
        Should be named forcePaint?
        Use me to force a repaint, rather than calling OnPaint directly.
        """
        self.OnPaint(None)
        return(True)

    #------------------------------------------------------------------
    # Internal painters
    #------------------------------------------------------------------

    def __realToLocal(self, x, y):
        """
        (Internal)
        Convert real genomic coords to local pixel coordinates.
        local to the location of the genome line.
        """
        return((self.w * ((x-self.lbp) / self.deltaf), (self.h - 30) + y))

    def __localToReal(self, sx, sy):
        """
        opposite of realToLocal()
        """
        x = 0
        y = 0
        return((x,y))

    def __getTextExtents(self, text):
        """
        (Internal)
        returns a bounding box of the text return = (x,y,w,h)
        """
        return(self.ctx.text_extents(text)[:4])

    def __drawTrackBackground(self, track_location, track_type):
        """
        track_location is the bottom edge of the track block
        """
        # get an available track slot
        self.__setPenColour( (0.95,0.95,0.95) )
        base_loc = self.__realToLocal(0, track_location)
        self.ctx.rectangle(0, base_loc[1]-opt.track.height_px[track_type], self.w, opt.track.height_px[track_type]-2) # 30 = half genomic track size
        self.ctx.fill()
        return( (0, base_loc[1]-opt.track.height_px[track_type], self.w, opt.track.height_px[track_type]-2) )

    def __drawTrackGraph(self, track_data, scaled=True, min_scaling=2):
        """
        **Arguments**
            track_data
                must be some kind of array/iterable, with a 1px resolution.

            scaled (True|False)
                scale the data vertically for the available track
                height?

            min_scaling
                only works if scaled = True,
                sets it so that a height of 1 is not expanded to the
                full height of the track. Instead the track will be scaled to
                this value as a minimum.
        """
        if scaled:
            track_max = max(track_data["array"])

            if min_scaling and track_max < min_scaling:
                scaling_value = min_scaling / float(opt.track.height_px["graph"])
            else:
                scaling_value = track_max / float(opt.track.height_px["graph"])
            new_array = track_data["array"] / scaling_value
        else:
            new_array = track_data["array"]

        colbox = self.__drawTrackBackground(track_data["track_location"], "graph")
        self.__setPenColour( (0,0,0) )
        self.ctx.set_line_width(0.5)
        coords = []
        lastpx = -1
        for index, value in enumerate(new_array):
            loc = self.__realToLocal(self.lbp + index, track_data["track_location"])
            #if int(loc[0]) > lastpx: # this means only draw one per pixel
            #    lastpx = loc[0]
            #    coords.append( (index, loc[1] - 30 - value)) # +30 locks it to the base of the track
            coords.append( (index, loc[1] - value)) # +30 locks it to the base of the track

        self.ctx.move_to(coords[0][0], coords[0][1]) # start x,y
        for item in coords:
            self.ctx.line_to(item[0], item[1])

        self.ctx.stroke()
        self.__drawText(0, loc[1] - 15 , opt.graphics.font, track_data["name"])

        return(colbox)# collision box dimensions

    def __drawTrackGraph_split_strand(self, track_data, scaled=True, min_scaling=100):
        """
        **Purpose**
            Similar to Graph, but draws the track into a top strand and a bottom strand.
            Although I could do this by generalising __drawTrackGraph() I feel it makes more
            sense not to re-use the code.

        **Arguments**
            track_data
                must be a dictionary of the form {"+": array(), "-": array()}
                where the arrays are iterables containing x bp resolution
                data.

            scaled (True|False)
                scale the data vertically for the available track
                height?

            min_scaling
                only works if scaled = True,
                sets it so that a height of 1 is not expanded to the
                full height of the track. Instead the track will be scaled to
                this value as a minimum.
        """
        assert "+" in track_data["array"], "__splitgraph data is missing + strand"
        assert "-" in track_data["array"], "__splitgraph data is missing - strand"

        half_way_point = opt.track.height_px["graph"] / 2

        if scaled:
            track_max = max([max(track_data["array"]["+"]), max(track_data["array"]["-"])])

            if min_scaling and track_max < min_scaling:
                scaling_value = min_scaling / float(half_way_point)
            else:
                scaling_value = track_max / float(half_way_point)
            # only works if numpy array?
            new_f_array = track_data["array"]["+"] / scaling_value # okay numpy can be sweet
            new_r_array = track_data["array"]["-"] / scaling_value # okay numpy can be sweet
        else:
            new_f_array = track_data["array"]["+"]
            new_r_array = track_data["array"]["-"]

        colbox = self.__drawTrackBackground(track_data["track_location"], "graph")

        for i, s in enumerate([new_f_array, new_r_array]):
            # + strand:
            if i == 0:
                self.__setPenColour( (0.8,0,0) )
            elif i == 1:
                self.__setPenColour( (0,0,0.8) )
            self.ctx.set_line_width(0.5)
            coords = []
            lastpx = -1
            for index, value in enumerate(s):
                loc = self.__realToLocal(self.lbp + index, track_data["track_location"])
                # get the middle:
                middle = loc[1] - half_way_point
                if i == 0:
                    coords.append( (index, middle + value)) # +30 locks it to the base of the track
                elif i == 1:
                    coords.append( (index, middle - value)) # +30 locks it to the base of the track

            self.ctx.move_to(coords[0][0], coords[0][1]) # start x,y
            for index, item in enumerate(coords):
                self.ctx.line_to(item[0], item[1])
            self.ctx.stroke()

        # - strand

        self.__drawText(0, loc[1] - 15 , opt.graphics.font, track_data["name"])
        return(colbox)

    def __drawTrackSpot(self, track_data, **kargs):
        """
        **Arguments**
            track_data
                a result from track.get_locations()

            colour
                a float colour (r, g, b, [a]), ranging 0..1
        """
        colbox = self.__drawTrackBackground(track_data["track_location"], "spot")

        sc = self.__realToLocal(0, track_data["track_location"])

        self.__setPenColour((1,1,1))
        self.ctx.rectangle(0, sc[1]-(opt.track.spot_pixel_radius*3), self.w, (opt.track.spot_pixel_radius*2)-2) # 30 = half genomic track size
        self.ctx.fill()

        colour = opt.track.spot_default_colour
        if "colour" in kargs:
            colour = kargs["colour"]
        self.__setPenColour(colour)

        for item in track_data["array"]:
            if opt.track.spot_shape == "circle":
                centre_point = (item["left"] + item["right"]) / 2
                sc = self.__realToLocal(centre_point, track_data["track_location"])
                self.ctx.arc(sc[0], sc[1]-(opt.track.spot_pixel_radius * 2), opt.track.spot_pixel_radius, 0, 2 * math.pi)
            elif opt.track.spot_shape == "triangle":
                pass

        if opt.track.spot_filled:
            self.ctx.fill()
        else:
            self.ctx.stroke()

        self.__drawText(0, sc[1]-17 , opt.graphics.font, track_data["name"])
        return(colbox)

    def __drawText(self, x, y, font, text, size=12, colour=(0,0,0), style=None):
        """
        (Internal - helper)
        Draw text to the screen, font is a string for the font name.
        Does not check for availability of the font.
        Will accept several styles, see constants for details.

        """
        self.__setPenColour(colour)
        self.ctx.select_font_face(font, txtToCairoA[style], txtToCairoB[style])
        self.ctx.move_to(x, y)
        self.ctx.show_text(text)

    def __drawGene(self, data):
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
        #print "t:", ((data["left"]-self.lbp) / self.deltaf), ((data["right"]-self.lbp) / self.deltaf)

        posLeft = self.__realToLocal(data["loc"]["left"], 0)[0]
        posRight = self.__realToLocal(data["loc"]["right"], 0)[0]
        self.ctx.set_line_width(1)
        self.__setPenColour(opt.graphics.gene_colour)

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
        coords.append(self.__realToLocal(data["loc"]["left"], 0))
        coords.append(self.__realToLocal(data["loc"]["left"], + opt.graphics.gene_height))
        for c in tc:
            if c["t"] == "cdss":
                current_offset = opt.graphics.cds_height
                coords.append(self.__realToLocal(c["c"], + opt.graphics.gene_height))
                coords.append(self.__realToLocal(c["c"], + opt.graphics.cds_height))
            elif c["t"] == "cdse":
                current_offset = opt.graphics.gene_height
                coords.append(self.__realToLocal(c["c"], + opt.graphics.cds_height))
                coords.append(self.__realToLocal(c["c"], + opt.graphics.gene_height))
            elif c["t"] == "es":
                coords.append(self.__realToLocal(c["c"], 0))
                coords.append(self.__realToLocal(c["c"], + current_offset))
            elif c["t"] == "ee":
                coords.append(self.__realToLocal(c["c"], + current_offset))
                coords.append(self.__realToLocal(c["c"], 0))

        coords.append(self.__realToLocal(data["loc"]["right"], + opt.graphics.gene_height))
        coords.append(self.__realToLocal(data["loc"]["right"], 0))

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
            loc = self.__realToLocal(data["loc"]["left"], 0)
            # arrow.
            self.ctx.move_to(loc[0], loc[1]-20)
            self.ctx.line_to(loc[0], loc[1]-20-opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0]+opt.graphics.arrow_width_px, loc[1]-20)
            self.ctx.line_to(loc[0], loc[1]-20+opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0], loc[1]-20)
            self.ctx.fill()
            self.__drawText(loc[0]+opt.graphics.arrow_width_px+3, loc[1]-20+opt.graphics.arrow_height_px, "Arial", data["name"], 12)
        elif data["strand"] == "-":
            loc = self.__realToLocal(data["loc"]["right"], 0)
            # arrow.
            self.ctx.move_to(loc[0], loc[1]+20)
            self.ctx.line_to(loc[0], loc[1]+20-opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0]-opt.graphics.arrow_width_px, loc[1]+20)
            self.ctx.line_to(loc[0], loc[1]+20+opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0], loc[1]+20)
            self.ctx.fill()
            self.__drawText(loc[0]+opt.graphics.arrow_width_px+3, loc[1]+20+opt.graphics.arrow_height_px, "Arial", data["name"], 12)
        else:
            raise ErrorInvalidGeneDefinition


        if opt.draw.single_midline_in_introns: # draw a single line through the gene
            # this looks best when the genome is not being drawn.
            leftmost = self.__realToLocal(data["loc"]["left"], 0)
            rightmost = self.__realToLocal(data["loc"]["right"], 0)
            self.__setPenColour(opt.graphics.gene_colour)
            self.ctx.set_line_width(1)
            self.ctx.move_to(leftmost[0], leftmost[1])
            self.ctx.line_to(rightmost[0], rightmost[1])
            self.ctx.stroke()

        if opt.draw.chevrons_inside_introns:
            pass

        if opt.draw.braces_between_exons:
            pass

        return(True)

    def __drawTrackBar(self, track_data, min_scale=1):
        """
        draw a 'bar' format track

        **Arguments**
            track_data
                must be some kind of array/iterable, with a nbp:1px resolution.

        """
        track_max = max(track_data["array"]) # bartrack must be normalised
        if track_max < min_scale:
            track_max = min_scale
        new_array = track_data["array"]

        posLeft = self.__realToLocal(self.lbp, track_data["track_location"])
        posRight = self.__realToLocal(self.rbp, track_data["track_location"])

        colbox = self.__drawTrackBackground(track_data["track_location"], "bar")

        self.ctx.set_line_width(10)

        currValue = new_array[0]
        col = 1.0 - (currValue / track_max)
        self.__setPenColour( (col,col,col) )
        self.ctx.move_to(posLeft[0], posLeft[1]-9) # start x,y

        for index, value in enumerate(new_array):
            fraction_along_array = index / len(new_array)

            if value != currValue:
                # move to the new position-1 and complete the line
                self.ctx.line_to(index, posLeft[1]-9)
                self.ctx.stroke()

                #change the colour to the new value:
                col = 1.0 - (value / track_max)
                self.__setPenColour( (col,col,col) )

                # move to the new start of the lien:
                self.ctx.move_to(index, posLeft[1]-9)

                # update the currValue
                currValue = value
                #print index, col, len(new_array)

        self.ctx.line_to(posRight[0], posRight[1]-9)
        self.ctx.stroke()

        self.__drawText(0, posRight[1]-17 , opt.graphics.font, track_data["name"])
        return(colbox)

    def __setPenColour(self, colour):
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
        this shouldn't be here?

        returns a collision type (see boundbox)
        or FALSE
        """
        for box in self.colBoxes:
            c = box.collideB((x,y))
            if c:
                print c["type"]
                return(c["type"])
