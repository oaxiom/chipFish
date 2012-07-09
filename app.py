
import sys, os 

import opt, gDraw

from glbase_wrapper import glload, track, location, format, flat_track
from error import *

class app():
    """
    Inheritence was split off a while ago.
    """
    def __init__(self, genome=None):
        """
        Any preamble set-up
        """
        self.current_genome = None

    def restart(self, tracklist=None):
        """
        see startup()
        
        Any clearing out can be performed here.
        This is not currently in use and is untested.
        This method would be used to rebind a tracklist
        """
        # Presumably here I should clean-out the tracklists...
        # Whoops!
        
        # reload through the normal route.
        self.startup(genome, tracklist)

    def __do_options(self, options, mode):
        """
        sort out the options. returns a dictionary of option:value pairs
        
        options should be a string containing option=value. Several options can be passed.
        mode will set default options for the track type
        """
        if not "=" in options:
            return({}) # No valid parseable options
        
        # Standardise input string:
        s = options.strip(" ").replace("\t", " ").replace("  ", " ").replace("\n", "").replace("\r", "")
        t = s.split(" ")
        res = {}
        for i in t:
            if "=" in i: # ignore mangled options
                opt, val = i.split("=")
                try: # try to cooerce ints and floats and bools
                    if "." in i.split("=")[1]: # try float
                        res[i.split("=")[0]] = float(i.split("=")[1])
                    elif val == "True":
                        res[opt] = True
                    elif val == "False":
                        res[opt] = False
                    else:
                        res[i.split("=")[0]] = int(i.split("=")[1])
                except ValueError:
                    res[i.split("=")[0]] = i.split("=")[1]
        
        return(res)

    def startup(self, tracklist=None):
        """
        startup or restart the server.

        Use this to get started proper and change genomes etc.

        **Arguments**
            tracklist
                A text file containing the path to the track files
                Not Implemented
        """
        self.draw = gDraw.gDraw()

        if tracklist: # Although it doesn't really make much sense not to supply a tracklist
            oh = open(tracklist, "rU")
            track_path = os.path.dirname(tracklist)
            mode = None
            for line in oh:
                if "#" in line:
                    # Options should go here.
                    pass # Anything useful to do with these?
                else:
                    if ":" in line: 
                        # The order of the following modes is important - for example, "bed:" will also
                        # collect "macs_bed" so "macs_bed" must go first.
                        # There can only be one mode
                        if "kde_track:" in line: # Must go before "track"
                            mode = "kde_track"
                        elif "track:" in line:
                            mode = "track"
                        elif "macs_bed:" in line:
                            mode = "macs_bed"
                        elif "bed:" in line:
                            mode = "bed"
                        elif "genome:" in line:
                            mode = "genome"
                        elif "flat:" in line:
                            mode = "flat"
                        else:
                            raise ErrorUnrecognisedTrackMode, mode
                            
                        # process options
                        options = self.__do_options(line.split(":")[-1], mode)
                    elif mode:
                        path = track_path
                        name = line.strip()
                        print options
                        if "abs_path" in options and options["abs_path"] == "True":
                            tail, head = os.path.split(name)
                            # The path will be relative to the path, not relative to chipFish. Which could be anywhere
                            # 
                            path = os.path.normpath(os.path.join(track_path, tail))
                            name = head
                            print path, name
                            
                        if name:
                            if mode == "track":
                                self.draw.bindTrack(track(filename=os.path.join(path, name)), options=options)
                                print "Bound Track:", os.path.join(path, name)
                            elif mode == "flat":
                                self.draw.bindTrack(flat_track(filename=os.path.join(path, name), bin_format="f"), track_type="graph", options=options)
                                print "Bound Track:", os.path.join(path, name)
                            elif mode == "kde_track":
                                self.draw.bindTrack(track(filename=os.path.join(path, name)), options=options, track_type="kde_graph")
                                print "Bound KDE Track:", os.path.join(path, name)
                            elif mode=="bed":
                                self.draw.bindTrack(genelist(filename=os.path.join(path, name), format=format.bed), options=options)
                                print "Bound Bed:", os.path.join(path, name) 
                            elif mode=="macs_bed":
                                f = format_minimal_bed
                                f["skiplines"] = 1 # Changes in recent version of macs bed.
                                self.draw.bindTrack(peaklist(filename=os.path.join(path, name), format=f), options=options)
                                print "Bound MACS bed file:", os.path.join(path, name)
                            elif mode=="genome": # must be a glb
                                self.draw.bindTrack(glload(os.path.join(path, name)))
                                print "Bound genome:", name 
            oh.close()

        self.draw.setViewPortSize(opt.draw.view_port_width)
        self.draw.setLocation(loc=location(loc="chr1:172724244-172859108")) # Load a dummy location, my favourite gene Stat3