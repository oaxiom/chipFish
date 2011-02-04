

import gDraw

from glbase_wrapper import glload, track, location

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
            for line in oh:
                if "#" in line:
                    pass # Anything useful to do with these?
                else:
                    self.draw.bindTrack(track(filename=os.pth.join(path, line)))
            oh.close()

        self.draw.setLocation(loc=location(loc="chr1:172724244-172859108")) # Load a dummy location.
