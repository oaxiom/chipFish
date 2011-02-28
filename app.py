
import sys, os 

import gDraw

from glbase_wrapper import glload, track, location, peaklist, format_bed
from error import *

class app():
    """
    This just inherits at the moment, but later I can use it to split off the gui and
    drawing elements, and tidy up gDraw to make a lot more public functions

    """
    def __init__(self, genome=None):
        """
        Any preamble set-up
        """
        self.current_genome = None

    def restart(self, genome="mm9", tracklist=None):
        """
        see startup()
        """
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

    def startup(self, genome="mm9", tracklist=None):
        """
        startup or restart the server.

        Use this to get started proper and change genomes etc.

        **Arguments**
            genome [mm8|mm9]
                load a supported genome.

            tracklist
                A text file containing the path to the track files
                Not Implemented
        """
        if genome == "mm8":
            self.g = glload("../Tracks/mm8_refGene.glb")
        elif genome == "mm9":
            self.g = glload("data/mm9_refGene.glb")
        else:
            raise ErrorGenomeNotSupported, genome
        self.current_genome = genome

        self.draw = gDraw.gDraw(self.g) # drawer must have access to the genome

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
                        if "kde_track" in line: # Must go before "track"
                            mode = "kde_track"
                        elif "track" in line:
                            mode = "track"
                        elif "bed" in line:
                            mode = "bed"
                        else:
                            raise ErrorUnrecognisedTrackMode, mode
                            
                        # process options
                        options = self.__do_options(line.split(":")[-1], mode)
                    elif mode:
                        name = line.replace("\n", "").replace("\r", "").replace("\t", "")
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

            oh.close()

        self.draw.setLocation(loc=location(loc="chr1:172724244-172859108")) # Load a dummy location.
