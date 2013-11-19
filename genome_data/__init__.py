"""

Sort out the available genomes

"""

import sys, os, glob
#sys.path.append(os.path.realpath("../")) # only required for testing
from glbase_wrapper import glload

class genomes:
    def __init__(self):
        """
        Work out available genome *.glb
        """
        self.datadir = os.path.dirname(os.path.realpath(__file__))
        self.genomes = []
        
        for genome in glob.glob(os.path.join(self.datadir, "*.glb")):
            self.genomes.append(os.path.split(genome)[1].split(".")[0])
        
    def get_genome(self, genome):
        assert genome in self.genomes, "genome '%s' is not available"
               
        return(glload(os.path.join(self.datadir, "%s.glb" % genome)))
        
    def __contains__(self, value):
        return(value in self.genomes)
        
if __name__ == "__main__":
    g = genomes()
    print g.genomes
    
    print "mm9" in g
    print "mm9_refGene" in g
    
    print g.get_genome("mm10_refGene")