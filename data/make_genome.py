"""
'make' a genome

turn refGene table (from mm8) into a genome.
"""

import sys, csv

sys.path.append("../")  # get the glbase_wrapper
import glbase_wrapper

# Build the mm8 genome.
mm8_refGene_format = {
    "strand": 3,
    "loc": {"code": "location(chr=column[2], left=column[4], right=column[5])"},
    "cds_loc": {"code": "location(chr=column[2], left=column[6], right=column[7])"},
    "exonCount": 7,
    "exonStarts": {"code": "[int(x) for x in column[9].strip(\",\").split(\",\")]"},
    "exonEnds": {"code": "[int(x) for x in column[10].strip(\",\").split(\",\")]"},
    "proteinID": 10,
    #"alignID": line[11],
    "name": 12,
    "refseq": 1,
    "dialect": csv.excel_tab
    } # this is how it will come out of the db.

mm8 = glbase_wrapper.genome(name="mm8", filename="mm8_refGene.tsv", format=mm8_refGene_format)
print mm8
mm8.save("mm8_refGene.glb")
print mm8._findByLabel("name", "Nanog")
