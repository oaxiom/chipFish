"""
--------------------------
chipFish, part of chipFish

(c) 2008-2018 oAxiom

--------------------------

Main app entry

NOTES:
------

"""

import sys, os, copy

import opt, gDraw


from glbase_wrapper import glload, track, location, format, flat_track, genelist, genome_sql
from error import *
from genome_data import genomes

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
                    res[i.split("=")[0]] = i.split("=")[1] # String literal

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
                if not line.strip(): # tolerate empty lines in spec sheet
                    continue

                if "#" in line[0]:
                    # Options should go here?
                    pass # Anything useful to do with these?
                else:
                    if ":" in line:
                        # The order of the following modes is important - for example, "bed:" will also
                        # collect "macs_bed" so "macs_bed" must go first.
                        # There can only be one mode
                        # Gods, just do strip() and regex it!?!?!!!!
                        if "kde_track:" in line: # Must go before "track"
                            mode = "kde_track"
                        elif "split_track" in line:
                            mode = "split_track"
                        elif "track:" in line:
                            mode = "track"
                        elif "macs_bed:" in line: # must go before bed
                            mode = "macs_bed"
                        elif "bed:" in line:
                            mode = "bed"
                        elif "genome:" in line:
                            mode = "genome"
                        elif "flat:" in line:
                            mode = "flat"
                        elif 'genome_sql:' in line:
                            mode = 'genome_sql'
                        else:
                            raise ErrorUnrecognisedTrackMode(mode)

                        # process options
                        options = self.__do_options(line.split(":")[-1], mode)
                    elif mode:
                        path = track_path
                        tt = line.strip().split()
                        name = tt[0]
                        track_opts = {}
                        if len(tt) > 1:
                            track_opts = self.__do_options(tt[1], mode)

                        if "abs_path" in options and options["abs_path"] == "True":
                            tail, head = os.path.split(name)
                            # The path will be relative to the path, not relative to chipFish. Which could be anywhere
                            #
                            path = os.path.normpath(os.path.join(track_path, tail))
                            name = head
                            #print path, name

                        if name:
                            options = copy.deepcopy(options)
                            options.update(track_opts)

                            if mode == "track":
                                self.draw.bindTrack(track(filename=os.path.join(path, name)), options=options)
                            elif mode == "split_track":
                                self.draw.bindTrack(track(filename=os.path.join(path, name)), options=options, track_type="graph_split_strand")
                            elif mode == "flat":
                                self.draw.bindTrack(flat_track(filename=os.path.join(path, name)), track_type="graph", options=options)
                            elif mode == "kde_track":
                                self.draw.bindTrack(track(filename=os.path.join(path, name)), options=options, track_type="kde_graph")
                            elif mode == "bed":
                                self.draw.bindTrack(genelist(filename=os.path.join(path, name), format=format.bed), options=options)
                            elif mode == "macs_bed":
                                f = format.minimal_bed
                                f["skiplines"] = 1 # Changes in recent version of macs bed.
                                self.draw.bindTrack(genelist(filename=os.path.join(path, name), format=f), options=options)
                            elif mode == "genome": # must be a glb
                                # First see if I can get it out of the pre-packaged genomes:
                                try:
                                    g = genomes()
                                    if name in g:
                                        self.draw.bindTrack(g.get_genome(name))
                                        continue
                                except AssertionError:
                                    # Okay, that did'nt work. see If I can get a file in this dir:
                                    self.draw.bindTrack(glload(os.path.join(path, name)))
                                    continue
                                # Okay, assume the user knows what they are doing and just grab the file they asked for:
                                self.draw.bindTrack(glload(name))
                            elif mode == "genome_sql":
                                self.draw.bindTrack(genome_sql(filename=os.path.join(path, name)))
            oh.close()

        self.draw.setViewPortSize(opt.draw.view_port_width)
        self.draw.setLocation(loc=location(loc="chr1:172724244-172859108")) # Load a dummy location, my favourite gene Stat3
