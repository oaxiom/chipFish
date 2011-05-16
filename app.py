
import sys, os 

import opt, gDraw

from glbase_wrapper import glload, track, location, peaklist, format_bed, format_minimal_bed
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
            try: # try to cooerce ints and floats
                if "." in i.split("=")[1]:
                    res[i.split("=")[0]] = float(i.split("=")[1])
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
            path = os.path.dirname(tracklist)
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
                        else:
                            raise ErrorUnrecognisedTrackMode, mode
                            
                        # process options
                        options = self.__do_options(line.split(":")[-1], mode)
                    elif mode:
                        name = line.replace("\n", "").replace("\r", "").replace("\t", "").replace(" ", "")
                        if name:
                            if mode=="track":
                                self.draw.bindTrack(track(filename=os.path.join(path, name)), options=options)
                                print "Bound Track:", os.path.join(path, name)
                            elif mode == "kde_track":
                                self.draw.bindTrack(track(filename=os.path.join(path, name)), options=options, track_type="kde_graph")
                                print "Bound KDE Track:", os.path.join(path, name)
                            elif mode=="bed":
                                self.draw.bindTrack(peaklist(filename=os.path.join(path, name), format=format_bed), options=options)
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