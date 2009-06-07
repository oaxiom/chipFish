"""
'make' a genome

build a refGene like list into a genome suitable for chipfish to display.
"""

import sys

# find glbase and import it.
sys.path.append("../../../glbase") # in the final version it would be in ".."

# modify the start-up of glbase
import config # some of the module names will clash...
config.REnv.libraries_to_test = []

# get the libraries I want.
from flags import *
from genome import genome

# Build the mm8 genome.
illuminaFormat = {"loc": 0, "symbol": 1, "refseq": 5, "entrez": 6, "strand": 7, "name": 8}

mm8 = genome("mm8", "illumina_remapped.csv", format_override=illuminaFormat)
print mm8
mm8.save("mm8.glb")
