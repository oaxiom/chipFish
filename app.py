

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

        if genome == "mm8":
            #self.draw.bindTrack(peaklist(filename="data/Matches_ESEsrrb_esrrf_5bp_soxf.bed", format=format_bed))
            #self.draw.bindTrack(peaklist(filename="../../Final Peak Lists/ESEsrrb_peaks.csv", format=format_chip_lists))
            #self.draw.bindTrack(peaklist(filename="../../Final Peak Lists/ESSox2_peaks.csv", format=format_chip_lists))
            #self.draw.bindTrack(peaklist(filename="../../Final Peak Lists/ESOct4_peaks.csv", format=format_chip_lists))
            #self.draw.bindTrack(track(filename="../Tracks/ES_H3K4me3.trk", name="ES H3K4me3"))
            #self.draw.bindTrack(track(filename="../Tracks/ES_H3K27me3.trk", name="H3K27me3"))
            #self.draw.bindTrack(track(filename="../Tracks/NS_CTCF.trk", name="ESCTCF"), track_type="graph_split_strand")
            #self.draw.bindTrack(track(filename="../Tracks/ESCTCF.trk", name="ESCTCF"))
            #self.draw.bindTrack(track(filename="../Tracks/ES_H3K36me3.trk", name="p300"))
            #self.draw.bindTrack(track(filename="data/TcMash1_new.trk", name="Telencephalon Mash1 ChIP-seq"))
            #self.draw.bindTrack(track(filename="data/NS_H3K4me3.trk", name="NS5 H3K4me3"))
            #self.draw.bindTrack(track(filename="data/NS_H3K36me3.trk", name="NS5 H3K36me3"), track_type="bar")
            #self.draw.bindTrack(track(filename="../Tracks/ES_H3K36me3.trk", name="H3K36me3"))
            #self.draw.bindTrack(track(filename="data/MEF_H3K4me3.trk", name="MEF H3K4me3"))

            self.draw.bindTrack(track(filename="../Tracks/TCD4_naive_H3k4me3.trk", name="TCD4 Naive H3K4me3"))
            self.draw.bindTrack(track(filename="../Tracks/TCD4_naive_H3k27me3.trk", name="TCD4 Naive H3K27me3"))
            #self.draw.bindTrack(track(filename="../Tracks/TCD4_Th1_H3k4me3.trk", name="TCD4 h1 H3K4me3"))
            self.draw.bindTrack(track(filename="../Tracks/TCD4_Th1_H3k27me3.trk", name="TCD4 h1 H3K27me3"))
            #self.draw.bindTrack(track(filename="../Tracks/TCD4_Th2_H3k4me3.trk", name="TCD4 h2 H3K4me3"))
            self.draw.bindTrack(track(filename="../Tracks/TCD4_Th2_H3k27me3.trk", name="TCD4 h2 H3K27me3"))
            #self.draw.bindTrack(track(filename="../Tracks/ESCTCF.trk", name="ESCTCF"))
        elif genome == "mm9":
            # mm9
            print track(filename="../Tracks/PEC_pIL10_aStat3.trk", name="PEC +IL10 anti-Stat3")

            self.draw.bindTrack(track(filename="../Tracks/PEC_pIL10_aStat3.trk", name="PEC +IL10 anti-Stat3"))
            #self.draw.bindTrack(track(filename="../Tracks/PEC_pIL10_Input.trk", name="PEC +IL10 Input"))
            self.draw.bindTrack(track(filename="../Tracks/PEC_Untre_aStat3.trk", name="PEC -IL10 anti-Stat3"))
            #self.draw.bindTrack(track(filename="../Tracks/PEC_Untre_Input.trk", name="PEC -IL10 Input"))

            # paired end RNA-seq tags
            self.draw.bindTrack(track(filename="../Tracks/PEC_Untre.trk", name="PEC -IL10"), track_type="graph_split_strand")
            self.draw.bindTrack(track(filename="../Tracks/PEC_pIL10.trk", name="PEC +IL10"), track_type="graph_split_strand")

        self.draw.setLocation(loc=location(loc="chr1:172724244-172859108")) # Load a dummy location.

