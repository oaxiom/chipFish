"""

This is a cli interface into chipfish.

I had to make this as the python bindings for wxwidgets failed to build
on the Mac, so I needed a way to output the images.

This is ripped out of gDraw and others.

"""

from __future__ import division

import sys, os, math, time

import opt

from error import *
from constants import *
from boundbox import bbox
from operator import itemgetter
from glbase_wrapper import location, track, utils
from data import *
from ruler import ruler
from log import log

MAX_TRACKS = 100 # maximum number of tracks This is kind of ludicrously high. 

import cairo, numpy

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
        self.paintQ = []
        self.tracks = [] # currently visible tracks
        self.trackBoxes = [False for n in xrange(MAX_TRACKS)] # list of TrackBoxes in use.

        self.chromosome = "1"
        self.lbp = 1
        self.rbp = 1000

        # set up dummy values for the view
        self.w = 100
        self.h = 200
        self.ctx = None
        self.last_guess_height = -1
        
        self.ruler = ruler(self.lbp, self.rbp, (0, self.w), "bp", True)

    def paint(self, cairo_context):
        """
        **Purpose**
            paint the current display on a cairo context
        
        **Arguments**
            cairo_context
                I need a valid cairo_context to draw to
        """
        self.ctx = cairo_context

        draw_modes_dict = {
            "gene": self.__drawGene,
            "lncRNA": self.__drawGene,
            "microRNA": self.__drawGene,
            "graph": self.__drawTrackGraph,
            "kde_graph": self.__drawTrackGraph,
            "graph_split_strand": self.__drawTrackGraph_split_strand,
            "bar": self.__drawTrackBar,
            "spot": self.__drawTrackSpot,
            "genome": self.__drawGenome,
            "repeats": self.__drawRepeats,
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

        draw_data = []
        # collect the data for the tracks and draw them on the screen.
        for track in self.tracks:
            #print track
            # basic data:
            item = {"type": track["type"],
                "track_location": track["track_location"],
                "name": track["data"]["name"], # This is problematic for genome() and genelist
                "options": track["options"]}

            # Big If statement coming:
            if track["type"] == "graph":
                item["array"] = track["data"].get_data("graph", location(loc=self.curr_loc),
                    resolution=self.bps_per_pixel, **track["options"])
                    
            elif track["type"] == "kde_graph":
                item["array"] = track["data"].get_data("graph", location(loc=self.curr_loc),
                    resolution=self.bps_per_pixel, kde_smooth=True, view_wid=self.w, **track["options"])
                    
            elif track["type"] == "bar":
                item["array"] = track["data"].get_data("bar", location(loc=self.curr_loc),
                    resolution=self.bps_per_pixel, **track["options"])
                    
            elif track["type"] == "spot":
                item["array"] = track["data"].get_data("spot", location(loc=self.curr_loc), **track["options"])
                
            elif track["type"] == "graph_split_strand":
                item["array"] = track["data"].get_data("graph", location(loc=self.curr_loc),
                    strand=True, resolution=self.bps_per_pixel, **track["options"])
                
            elif track["type"] == "genome":
                item["array"] = track["data"].get_data("genome", location(loc=self.curr_loc))
                #print item["array"]

            elif track["type"] == "repeats":
                item["array"] = track["data"].get_data(location(loc=self.curr_loc))

            draw_data.append(item)

        # If we need to lock the tracks:
        if opt.track.lock_scales:
            # Fid the max_values for all of the tracks and then load min_scaling 
            scale = opt.track.min_scale
            for item in draw_data:
                mmax = 0
                if item["type"] in ["graph", "kde_graph"]:
                    mmax = max(item["array"])
                elif item["type"] == "graph_split_strand": # got two max's to do.
                    mmax = max([max(item["array"]["+"]), max(item["array"]["-"])])
                
                if mmax > scale: # see if we need to adjust scale
                    scale = mmax
                    
            # Now set all tracks to the same scale.
            for item in draw_data:
                if item["type"] in ["graph", "kde_graph", "graph_split_strand"]:
                    item["options"]["min_scaling"] = scale
                                
        # And finally draw:
        for item in draw_data:
            #print item
            colbox = draw_modes_dict[item["type"]](item, **item["options"])
             
            # the collision boxes are not used. But I suppose in future...
            if colbox:
                self.__col_boxs.append(bbox(colbox, track, "track"))

        #self.__col_boxs.append(bbox(self.__drawRuler(), None, "ruler")) # Unused at the moment, but may be useful later.
        if opt.ruler.draw:
            self.ruler.draw(self.ctx, (0,0))
            
        if opt.draw.genomic_location:
            self.__drawText(5, opt.ruler.height_px + 22, opt.graphics.font, str(self.curr_loc), size=opt.draw.genomic_location_font_size)
            
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

    def bindTrack(self, track, options=None, track_type=None):
        """
        bind a drawing track extra to genome.
        
        valid track types:
        %s
        """ % valid_track_draw_types
        # if no track_type try to guess from the track object
        if options and 'track_type' in options:
            track_type = options['track_type']
            
        if not track_type:
            track_type = track._default_draw_type

        if track_type not in valid_track_draw_types:
            raise ErrorTrackDrawTypeNotFound, track_type

        if not options:
            options = {} # solves a lot of problems if it is an empty dict

        self.tracks.append({"data": track, "track_location": 
            self.__getNextTrackBox(track_type), "type": track_type,
            "options": options})

    def __getNextTrackBox(self, track_type):
        """
        get the next available track location and return the base (y-axis) of the block.
        """
        currentLoc = 0
        for index, track in enumerate(self.trackBoxes):
            if track:
                currentLoc += opt.track.height_px[track]
            if not track:
                self.trackBoxes[index] = track_type
                return(-(index + (currentLoc))) 

    def setViewPortSize(self, width):
        """
        **Purpose**
            set the size of the viewport
        
        **Arguments**
            width
                sets the width of the display
        """
        self.fullw = width # for blanking screen
        self.w = width - opt.graphics.right_border_width # a small right most border for editing trakcs, never implemented
        self.halfw = self.w / 2
        self.ruler.set(0, 100, (0, self.w))

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
            paint() or forceRedraw() correct.
        """
        if chromosome: # old-style assignation
            self.chromosome = str(chromosome)
            self.lbp = leftBasePair
            self.rbp = rightBasePair
            self.loc = location(chr=chromosome, left=self.lbp, right=self.rbp)
        elif loc: # new-style <location> assignation.
            self.chromosome = loc["chr"]
            self.lbp = loc["left"]
            self.rbp = loc["right"]
            self.loc = loc
        
        # sanity checking? Neccesary? yes
        # There's a bug here is the d[lbp, rbp] is less than w
        self.ruler.set(self.lbp, self.rbp, (0, self.w))

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

        scales_and_labels = {
            # 1 bp causes chipFish to crash :(
            10: "10bp",
            100: "100bp",
            1000: "1kbp",
            10000: "10kbp",
            100000: "100kbp",
            1000000: "1Mbp",
            10000000: "10Mbp",
            100000000: "100Mbp", # Human chr 1 is 250Mbp.
            1000000000: "1000Mbp", # Just in case there are some wierd genomes.
            10000000000: "1Gbp" # And like, some aliens or something 
            }

        for i, v in enumerate(scales_and_labels):
            if minim > abs(percent - v):
                minim = abs(percent - v)
                best = v

        posLeft = self.__realToLocal(self.rbp - best, -40)
        posRight = self.__realToLocal(self.rbp, -40)

        self.ctx.set_line_width(4.0)
        self.__setPenColour((0,0,0))
        
        self.ctx.move_to(posLeft[0]-20, opt.ruler.height_px + 8) # move 20px arbitrarily left
        self.ctx.line_to(posRight[0]-20, opt.ruler.height_px + 8)
        self.ctx.stroke()

        self.ctx.move_to(posLeft[0]-20, opt.ruler.height_px + 8 -4) # move 20px arbitrarily left
        self.ctx.line_to(posLeft[0]-20, opt.ruler.height_px + 8 +4)
        self.ctx.stroke()

        self.ctx.move_to(posRight[0]-20, opt.ruler.height_px + 8 -4) # move 20px arbitrarily left
        self.ctx.line_to(posRight[0]-20, opt.ruler.height_px + 8 +4)
        self.ctx.stroke()

        self.__drawText(posLeft[0], opt.ruler.height_px + 21, opt.graphics.font, scales_and_labels[best], size=opt.draw.scale_bar_fontsize)

    def exportImage(self, filename, scale=1, type=None):
        """
        **Purpose**

            export the current Cairo image as a 'type'
            the image type will be guessed from the filename.
            If no obvious extension is given then the string value in 'types'
            will be used.
            Finally, if that doens't make sense then it will default to a png
            
            <This is not implented - it only saves as a png>

        **Arguments**

            filename (Required)
                filename and path (and extension) to the file to save

            type (Optional)
                the type of file to save as, this will override the
                extension in the filename. If not type is given and
                the extension doesn't make any sense then a png will be
                used.
                
            scale (Optional, default=1)
                rescale the svg

        **Result**
            returns the actual filename used to save.
            and a file saved in filename
        """
        # guess the maximum height required:
        guess_height = abs(self.tracks[-1]["track_location"]) + opt.track.height_px[self.tracks[-1]["type"]] + opt.ruler.height_px + 30 + 32 # 32 is the 'chromosome %s' padding, 30 is some other padding I'm not really certain where it comes from...

        if not self.ctx or guess_height != self.last_guess_height or type in ("svg", "pdf"):
            # If the surface size changed between the last call and this, I need a newly sized surface.
            if type == "svg":
                self.surface = cairo.SVGSurface(filename, int(self.w*scale), int(guess_height*scale))
            elif type == "pdf": # not tested
                self.surface = cairo.PDFSurface(filename, int(self.w*scale), int(guess_height*scale))
            else: # get a png
                self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(self.w*scale), int(guess_height*scale))
            self.ctx = cairo.Context(self.surface)
            self.last_guess_height = guess_height

        self.h = guess_height

        # forceRedraw onto my surface.
        if scale != 1:
            self.ctx.scale(scale, scale) 
            
        self.paint(self.ctx)
        
        if scale != 1:
            self.ctx.scale(1.0/scale, 1.0/scale) 


        if type in ("svg", "pdf", "ps", "eps"):
            actual_filename = filename # suppose to mod, based on type
            self.surface.finish()
        else: # png
            # save image
            actual_filename = filename
            filehandle = open(actual_filename, "wb")
            self.surface.write_to_png(filehandle)
            filehandle.close()

        return(actual_filename)

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
        self.ctx.rectangle(0, base_loc[1]-opt.track.height_px[track_type]+4, self.w, opt.track.height_px[track_type]-4) # 30 = half genomic track size
        self.ctx.fill()
        return( (0, base_loc[1]-opt.track.height_px[track_type], self.w, opt.track.height_px[track_type]-2) )

    def __drawTrackGraph(self, track_data, scaled=True, min_scaling=opt.track.min_scale, clamp=1, 
        no_scaling=False, colour=None, name=None, **kargs):
        """
        **Arguments**
            track_data
                must be some kind of array, with a 1px resolution.

            scaled (True|False)
                scale the data vertically for the available track
                height?

            min_scaling
                only works if scaled = True,
                sets it so that a height of 1 is not expanded to the
                full height of the track. Instead the track will be scaled to
                this value as a minimum.
            
            clamp (default=1)
                clamp the display scale from <clamp> .. n
                rather than 0 .. n
                
            name (Optional, default=None)
                By default I will use the name of the genelist.
                If you want to rename the track then set options name="Name of track" 
        """

        data = track_data["array"]
        data = numpy.array(data, dtype=numpy.float32)

        if clamp:
            for i, v in enumerate(data):
                if v > clamp-1.0:
                    data[i] = v-clamp
                else:
                    data[i] = 0
        else:
            clamp = 0 # I need a clamp int later even if clamp=None

        track_max = max(data)
        track_min = min(data)

        if scaled:
            scaling_value = max(min_scaling+clamp, track_max) / float(opt.track.height_px["graph"])
            data = data / scaling_value
        #print ":", track_max, scaling_value, min_scaling,  max(min_scaling, track_max) 

        if opt.track.background:
            colbox = self.__drawTrackBackground(track_data["track_location"], "graph")
        else:
            colbox = []
        
        if not colour:
            self.__setPenColour( (0,0,0) )
        else:
            if isinstance(colour, str):
                self.__setPenColour(colour_lookup_name[colour])
            else:
                self.__setPenColour(colour) # this will probably not work.
        
        self.ctx.set_line_width(1.0)
        coords = []
        lastpx = -1
        for index, value in enumerate(data):
            loc = self.__realToLocal(self.lbp + index, track_data["track_location"])
            coords.append( (index, loc[1] - value)) # +30 locks it to the base of the track

        self.ctx.move_to(0, coords[0][1]) # start x,y
        
        for item in coords:
            self.ctx.line_to(item[0], item[1])
            
        if opt.track.filled:
            if clamp:
                loc = self.__realToLocal(1, track_data["track_location"])
            else:
                loc = self.__realToLocal(0, track_data["track_location"])
            #if clamp:
            self.ctx.line_to(item[0], loc[1] - min(data)) # move to the base line on the far right 
            self.ctx.line_to(0, loc[1] - min(data)) # the nth far left
            self.ctx.fill()
        else:
            self.ctx.stroke()
            
        if opt.track.draw_names:
            if not name:
                name = track_data["name"]
            self.__drawText(opt.track.label_fontsize, loc[1] - opt.track.height_px["graph"] + (opt.track.label_fontsize*2), 
                opt.graphics.font, name, size=opt.track.label_fontsize)
            
        if opt.track.draw_scales:
            self.__drawText(self.w - 10, loc[1] - 5, opt.graphics.font, 
                int(clamp), # Must be clamp
                size=opt.track.scale_bar_font_size, align="right", colour=(0,0,0))
            self.__drawText(self.w - 10, loc[1] - opt.track.height_px["graph"] + opt.track.scale_bar_font_size + 5, 
                opt.graphics.font, 
                int(max(min_scaling+clamp, track_max+clamp)), # You must add the clamp, otherwise you end up with maximums less than the clamp
                size=opt.track.scale_bar_font_size, align="right", colour=(0,0,0))

        return(colbox)# collision box dimensions

    def __drawTrackGraph_split_strand(self, track_data, scaled=True, min_scaling=opt.track.min_scale, 
        name=None, clamp=True, **kargs):
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

        track_max = max([track_data["array"]["+"].max(), track_data["array"]["-"].max()])
        track_min = min([track_data["array"]["+"].min(), track_data["array"]["-"].min()])

        if scaled:
            track_max = max([min_scaling, max(track_data["array"]["+"]), max(track_data["array"]["-"])])

            if min_scaling and track_max < min_scaling:
                scaling_value = min_scaling / float(half_way_point)
            else:
                scaling_value = track_max / float(half_way_point)
                
            # only works if numpy array?
            new_f_array = track_data["array"]["+"] / scaling_value # okay numpy can be sweet
            new_r_array = track_data["array"]["-"] / scaling_value 
        else:
            new_f_array = track_data["array"]["+"]
            new_r_array = track_data["array"]["-"]

        #print new_f_array
        #print track_data["array"]

        colbox = self.__drawTrackBackground(track_data["track_location"], "graph")

        for i, s in enumerate([new_f_array, new_r_array]): 
            if i == 0:# + strand:
                self.__setPenColour( (0, 0, 0.8) )
            elif i == 1:
                self.__setPenColour( (0.8, 0, 0) )
                
            self.ctx.set_line_width(2)
            
            # work out the list of screen-space coords
            #coords = []
            #lastpx = -1
            #for index, value in enumerate(data):
            #    loc = self.__realToLocal(self.lbp + index, track_data["track_location"])
            #    coords.append( (index, loc[1] - value)) # +30 locks it to the base of the track
            
            coords = []
            lastpx = -1
            for index, value in enumerate(s):
                loc = self.__realToLocal(self.lbp + index, track_data["track_location"])
                # get the middle:
                middle = loc[1] - half_way_point
                if i == 0:
                    coords.append( (index, middle - value)) # +30 locks it to the base of the track
                elif i == 1:
                    coords.append( (index, middle + value)) # +30 locks it to the base of the track

            self.ctx.move_to(coords[0][0], coords[0][1]) # start x,y
            for index, item in enumerate(coords):
                self.ctx.line_to(item[0], item[1]) # move along each coord
                
            if opt.track.filled:
                if clamp:
                    loc = self.__realToLocal(1, track_data["track_location"])
                else:
                    loc = self.__realToLocal(0, track_data["track_location"])

                self.ctx.line_to(item[0], middle) # move to the base line on the far right 
                self.ctx.line_to(0, middle) # the nth far left
                self.ctx.fill()
            else:
                self.ctx.stroke()

        if opt.track.draw_names:
            if not name:
                name = track_data["name"]
            self.__drawText(opt.track.label_fontsize, loc[1] - opt.track.height_px["graph"] + (opt.track.label_fontsize*2), 
                opt.graphics.font, name, size=opt.track.label_fontsize)
            
        if opt.track.draw_scales: # only track_max is drawn on split_graphs
            self.__drawText(self.w - 10, loc[1] - opt.track.height_px["graph"] + opt.track.scale_bar_font_size + 5, 
                opt.graphics.font, 
                int(max(min_scaling, track_max)), 
                size=opt.track.scale_bar_font_size, align="right", colour=(0,0,0))

        return(colbox)

    def __drawTrackSpot(self, track_data, **kargs):
        """
        **Arguments**
            track_data
                a result from track.get_locations()

            colour
                a float colour (r, g, b, [a]), ranging 0..1
        """
        #colbox = self.__drawTrackBackground(track_data["track_location"], "spot")

        sc = self.__realToLocal(0, track_data["track_location"])

        #self.__setPenColour((1,1,1))
        #self.ctx.rectangle(0, sc[1]-(opt.track.spot_pixel_radius*3), self.w, (opt.track.spot_pixel_radius*2)-2) # 30 = half genomic track size
        #self.ctx.fill()

        colour = opt.track.spot_default_colour
        if "colour" in kargs:
            colour = kargs["colour"]
        self.__setPenColour(colour)

        for item in track_data["array"]:
            if opt.track.spot_shape == "circle":
                centre_point = (item["left"] + item["right"]) / 2
                sc = self.__realToLocal(centre_point, track_data["track_location"])
                self.ctx.arc(sc[0], sc[1]-(opt.track.spot_pixel_radius * 2)-3, opt.track.spot_pixel_radius, 0, 2 * math.pi)
            elif opt.track.spot_shape == "triangle":
                pass

        if opt.track.spot_filled:
            self.ctx.fill()
        else:
            self.ctx.stroke()

        #self.__drawText(0, sc[1]-17 , opt.graphics.font, track_data["name"])
        return(None)

    def __drawTrackBar(self, track_data, min_scale=1, **kargs):
        """
        draw a 'bar-heatmap' format track

        **Arguments**
            track_data
                must be some kind of array/iterable, with a nbp:1px resolution.

        """
        #colbox = self.__drawTrackBackground(track_data["track_location"], "spot")

        sc = self.__realToLocal(0, track_data["track_location"])

        self.__setPenColour((1,1,1))
        #self.ctx.rectangle(0, sc[1], self.w, (opt.track.bar_height-1)) # 30 = half genomic track size
        #self.ctx.fill()

        colour = opt.track.bar_default_colour
        if "colour" in kargs:
            colour = kargs["colour"]
        self.ctx.set_line_width(0.5)
        self.__setPenColour(colour)

        for bar in track_data["array"]:           
            left = self.__realToLocal(bar['left'], track_data["track_location"])
            rite = self.__realToLocal(bar['right'], track_data["track_location"])
            
            self.ctx.move_to(left[0], sc[1]-opt.track.bar_height) 
            self.ctx.line_to(rite[0], sc[1]-opt.track.bar_height) 
            self.ctx.line_to(rite[0], sc[1]) 
            self.ctx.line_to(left[0], sc[1]) 
            self.ctx.fill()

        #self.__drawText(0, sc[1]-17 , opt.graphics.font, track_data["name"])
        return(None)

    def __drawTrackHeatmap(self, track_data, min_scale=1, **kargs):
        """
        draw a 'bar-heatmap' format track

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
        col = 1.0# - (currValue / track_max)
        self.__setPenColour( (col,col,col) )
        self.ctx.move_to(posLeft[0], posLeft[1]-9) # start x,y

        print new_array

        for index, value in enumerate(new_array):
            fraction_along_array = index / len(new_array)

            if value != currValue:
                # move to the new position-1 and complete the line
                self.ctx.line_to(index, posLeft[1]-9)
                self.ctx.stroke()

                #change the colour to the new value:
                col = 1.0# - (value / track_max)
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

    def __drawText(self, x, y, font, text, size=opt.graphics.default_font_size, colour=(0,0,0), style=None, align="left"):
        """
        (Internal - helper)
        Draw text to the screen, font is a string for the font name.
        Does not check for availability of the font.
        Will accept several styles, see constants for details.

        """
        self.__setPenColour(colour)
        self.ctx.select_font_face(font, txtToCairoA[style], txtToCairoB[style])
        self.ctx.set_font_size(size)
        if align == "left":
            self.ctx.move_to(x, y)
        elif align == "right":
            t = self.ctx.text_extents(str(text))
            self.ctx.move_to(x - t[2], y) 
        self.ctx.show_text(str(text))

    def __drawGenome(self, track_data, **kargs):
        # Draw a genome
        draw_modes_dict = { # In case I ever add other drawing modes...
            "gene": self.__drawGene,
            "lncRNA": self.__drawGene,
            "microRNA": self.__drawGene,
            }
    
        if opt.draw.double_lines_for_genome:
            self.__drawChr(None)

        # draw the genome:
        for item in track_data["array"]:
            draw_modes_dict[item["type"]](item, track_data)

    def __drawGene(self, data, track_data=None):
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
        if not track_data: # See if we are bound into a track slot.
            track_slot_base = 0
        else:
            track_slot_base = track_data["track_location"] - (opt.track.height_px['genome'] / 2.0)

        posLeft = self.__realToLocal(data["loc"]["left"], track_slot_base)[0]
        posBase = self.__realToLocal(data["loc"]["left"], track_slot_base)[1]
        posRight = self.__realToLocal(data["loc"]["right"], track_slot_base)[0]
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
        coords.append(self.__realToLocal(data["loc"]["left"], track_slot_base))
        coords.append(self.__realToLocal(data["loc"]["left"], track_slot_base + opt.graphics.gene_height))
        for c in tc:
            if c["t"] == "cdss":
                current_offset = opt.graphics.cds_height
                coords.append(self.__realToLocal(c["c"], track_slot_base + opt.graphics.gene_height))
                coords.append(self.__realToLocal(c["c"], track_slot_base + opt.graphics.cds_height))
            elif c["t"] == "cdse":
                current_offset = opt.graphics.gene_height
                coords.append(self.__realToLocal(c["c"], track_slot_base + opt.graphics.cds_height))
                coords.append(self.__realToLocal(c["c"], track_slot_base + opt.graphics.gene_height))
            elif c["t"] == "es":
                coords.append(self.__realToLocal(c["c"], track_slot_base))
                coords.append(self.__realToLocal(c["c"], track_slot_base + current_offset))
            elif c["t"] == "ee":
                coords.append(self.__realToLocal(c["c"], track_slot_base + current_offset))
                coords.append(self.__realToLocal(c["c"], track_slot_base))

        coords.append(self.__realToLocal(data["loc"]["right"], track_slot_base + opt.graphics.gene_height))
        coords.append(self.__realToLocal(data["loc"]["right"], track_slot_base))

        self.ctx.move_to(coords[0][0], coords[0][1])
        for index, item in enumerate(coords):
            self.ctx.line_to(item[0], item[1])
        self.ctx.move_to(coords[0][0], coords[0][1])

        coords.reverse()
        for item in coords:
            self.ctx.line_to(item[0], posBase + (posBase - item[1]))
            
        self.ctx.line_to(coords[0][0], coords[0][1]) # Finish so fill works 
        #self.ctx.stroke()
        self.ctx.fill()

        self.ctx.set_line_width(1)

        #---------------------------------------------------------------
        # Draw gene arrow
        '''
        if data["strand"] == "+": # top strand
            loc = self.__realToLocal(data["loc"]["left"], track_slot_base)
            # arrow.
            self.ctx.move_to(loc[0], loc[1]-20)
            self.ctx.line_to(loc[0], loc[1]-20-opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0] + opt.graphics.arrow_width_px, loc[1]-20)
            self.ctx.line_to(loc[0], loc[1]-20+opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0], loc[1]-20)
            self.ctx.fill()
            self.__drawText(loc[0]+opt.graphics.arrow_width_px+3, loc[1]-20+opt.graphics.arrow_height_px, "Arial", data["name"], size=opt.gene.font_size, style=opt.gene.font_style)
        elif data["strand"] == "-":
            loc = self.__realToLocal(data["loc"]["right"], track_slot_base)
            # arrow.
            self.ctx.move_to(loc[0], loc[1]+20)
            self.ctx.line_to(loc[0], loc[1]+20-opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0]-opt.graphics.arrow_width_px, loc[1]+20)
            self.ctx.line_to(loc[0], loc[1]+20+opt.graphics.arrow_height_px)
            self.ctx.line_to(loc[0], loc[1]+20)
            self.ctx.fill()
            self.__drawText(loc[0]+opt.graphics.arrow_width_px-13, loc[1]+22+opt.graphics.arrow_height_px, "Arial", data["name"], size=opt.gene.font_size, align="right", style=opt.gene.font_style)
        else:
            raise ErrorInvalidGeneDefinition
        '''

        if opt.draw.single_midline_in_introns: # draw a single line through the gene
            # this looks best when the genome is not being drawn.
            leftmost = self.__realToLocal(data["loc"]["left"], track_slot_base)
            rightmost = self.__realToLocal(data["loc"]["right"], track_slot_base)
            self.__setPenColour(opt.graphics.gene_colour)
            self.ctx.set_line_width(3)
            self.ctx.move_to(leftmost[0], leftmost[1])
            self.ctx.line_to(rightmost[0], rightmost[1])
            self.ctx.stroke()

        if opt.draw.chevrons_inside_introns:
            pass

        if opt.draw.braces_between_exons:
            pass

        return(True)

    def __drawRepeats(self, track_data, **kargs):
        # Draw a genome
        draw_modes_dict = { # In case I ever add other drawing modes...
            'LINE': self.__drawRepeat, 
            'LTR': self.__drawRepeat, 
            'SINE': self.__drawRepeat
            }

        self.__drawTrackBackground(track_data["track_location"], "repeats")
        
        one_third = ((opt.track.height_px['repeats']) / 4.0)
        # Three genome lines for the repeats
        self.ctx.set_line_width(1)
        self.__setPenColour((0.6, 0.6, 0.6))
        for i in (1.0, 2.0, 3.0):
            track_slot_base = track_data["track_location"] - (one_third * i)
            loc = self.__realToLocal(0, track_slot_base)
            self.ctx.move_to(0, loc[1])
            self.ctx.line_to(self.w, loc[1])
            self.ctx.stroke()

        # draw the genome:
        for item in track_data["array"]:
            # Should assert here if not in draw_modes
            if item['type'] not in draw_modes_dict:
                log.warning('type=%s not found in draw modes, skipping.')
                continue
            draw_modes_dict[item["type"]](item, track_data)

    def __drawRepeat(self, data, track_data=None):
        """
        draw Features of the type "Repeat"
        should be an dict containing the following keys:
        type: repeat
        loc: location span of gene
        strand: strand of gene
        Class: 
        """        
        one_third = ((opt.track.height_px['repeats']) / 4.0)
        # order = LINE, SINE, LTR
        if data['type'] == 'LINE':
            track_slot_base = track_data["track_location"] - one_third
        elif data['type'] == 'SINE':
            track_slot_base = track_data["track_location"] - (one_third * 2.0)
        elif data['type'] == 'LTR':
            track_slot_base = track_data["track_location"] - (one_third * 3.0)

        posLeft = self.__realToLocal(data["loc"]["left"], track_slot_base)[0]
        posBase = self.__realToLocal(data["loc"]["left"], track_slot_base)[1]
        posRight = self.__realToLocal(data["loc"]["right"], track_slot_base)[0]

        #---------------------------------------------------------------
        # draw repeat blocks.
        self.ctx.set_line_width(1)
        self.__setPenColour(opt.graphics.repeat_cols[data['type']])
        current_offset = opt.graphics.gene_height
        coords = []
        if data['strand'] == '+':
            coords.append(self.__realToLocal(data["loc"]["left"], track_slot_base))
            coords.append(self.__realToLocal(data["loc"]["left"], track_slot_base - opt.graphics.repeat_height))
            coords.append(self.__realToLocal(data["loc"]["right"], track_slot_base - opt.graphics.repeat_height))
            coords.append(self.__realToLocal(data["loc"]["right"], track_slot_base))
        elif data['strand'] == '-':
            coords.append(self.__realToLocal(data["loc"]["left"], track_slot_base))
            coords.append(self.__realToLocal(data["loc"]["left"], track_slot_base + opt.graphics.repeat_height))
            coords.append(self.__realToLocal(data["loc"]["right"], track_slot_base + opt.graphics.repeat_height))
            coords.append(self.__realToLocal(data["loc"]["right"], track_slot_base))

        self.ctx.move_to(coords[0][0], coords[0][1])
        for index, item in enumerate(coords):
            self.ctx.line_to(item[0], item[1])
        self.ctx.move_to(coords[0][0], coords[0][1])
            
        self.ctx.line_to(coords[0][0], coords[0][1]) # Finish so fill works 
        self.ctx.fill()

        self.ctx.set_line_width(1)
        self.__setPenColour(opt.graphics.repeat_cols[data['type']])

        #---------------------------------------------------------------
        # Draw label (strand is indicated by position relative to the genome line)

        if data["strand"] == "+": # top strand
            loc = self.__realToLocal(data["loc"]["left"], track_slot_base)
            # arrow.
            #self.ctx.move_to(loc[0], loc[1]-8)
            #self.ctx.line_to(loc[0], loc[1]-8-opt.graphics.repeat_arrow_height_px)
            #self.ctx.line_to(loc[0] + opt.graphics.arrow_width_px, loc[1]-8)
            #self.ctx.line_to(loc[0], loc[1]-8+opt.graphics.repeat_arrow_height_px)
            #self.ctx.line_to(loc[0], loc[1]-8)
            #self.ctx.fill()
            #self.__drawText(loc[0]+opt.graphics.arrow_width_px+3, loc[1]-8+opt.graphics.arrow_height_px, "Helvetica", data["name"], size=opt.graphics.repeat_label_font_size, style=opt.graphics.repeat_label_font_style)
            self.__drawText(loc[0], loc[1]-opt.graphics.repeat_height-2, "Helvetica", data["name"], size=opt.graphics.repeat_label_font_size, style=opt.graphics.repeat_label_font_style)
        elif data["strand"] == "-":
            loc = self.__realToLocal(data["loc"]["right"], track_slot_base)
            # arrow.
            #self.ctx.move_to(loc[0], loc[1]+8)
            #self.ctx.line_to(loc[0], loc[1]+8-opt.graphics.repeat_arrow_height_px)
            #self.ctx.line_to(loc[0]-opt.graphics.arrow_width_px, loc[1]+8)
            #self.ctx.line_to(loc[0], loc[1]+8+opt.graphics.repeat_arrow_height_px)
            #self.ctx.line_to(loc[0], loc[1]+8)
            #self.ctx.fill()
            #self.__drawText(loc[0]+opt.graphics.arrow_width_px-13, loc[1]+10+opt.graphics.arrow_height_px, "Helvetica", data["name"], size=opt.graphics.repeat_label_font_size, align="right", style=opt.graphics.repeat_label_font_style)
            self.__drawText(loc[0], loc[1]+opt.graphics.repeat_height+opt.graphics.repeat_label_font_size-2, "Helvetica", data["name"], size=opt.graphics.repeat_label_font_size, align="right", style=opt.graphics.repeat_label_font_style)
        else:
            raise ErrorInvalidGeneDefinition

        return(True)

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
