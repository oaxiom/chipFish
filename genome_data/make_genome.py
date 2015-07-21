"""
'make' a genome

turn refGene table into a genome.
"""

import sys, csv, os

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
    } # This is now in glbase >0.161.hg

if os.path.exists("mm8_refGene.tsv"):
    mm8 = glbase_wrapper.genome(name="mm8", filename="mm8_refGene.tsv", format=mm8_refGene_format)
    print mm8
    mm8.save("mm8_refGene.glb")
    print mm8._findByLabel("name", "Nanog")

format_mm9_refgene = {
        "loc": {"code": "location(chr=column[2], left=int(column[4]), right=int(column[5]))"},
        "strand": 3,
        "name": 12,
        "refseq": 1,
        #"tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"}, # Broken in this version of glbase
        "cds_loc": {"code": "location(chr=column[2], left=column[6], right=column[7])"},
        "exons_count": 8,
        "exonStarts": {"code": "[int(x) for x in column[9].strip(\",\").split(\",\")]"},
        "exonEnds": {"code": "[int(x) for x in column[10].strip(\",\").split(\",\")]"},
        "dialect" : csv.excel_tab
        } # This is now in glbase >0.161.hg

if os.path.exists("mm9_refGene.tsv"):
    mm9 = glbase_wrapper.genome(name="mm9", filename="mm9_refGene.tsv", format=format_mm9_refgene)
    print mm9
    mm9.save("mm9_refGene.glb")
    print mm9._findByLabel("name", "Nanog")

if os.path.exists("hg19_refGene.tsv"):
    hg19 = glbase_wrapper.genome(name="hg19", filename="hg19_refGene.tsv", format=format_mm9_refgene)
    print hg19
    hg19.save("hg19_refGene.glb")
    print hg19._findByLabel("name", "Nanog")

if os.path.exists("mm10_refGene.tsv"):
    mm10 = glbase_wrapper.genome(name="mm10", filename="mm10_refGene.tsv", format=format_mm9_refgene)
    print mm10
    mm10.save("mm10_refGene.glb")
    print mm10._findByLabel("name", "Nanog")
    
if os.path.exists("hg38_refGene.tsv"):
    mm10 = glbase_wrapper.genome(name="mm10", filename="hg38_refGene.tsv", format=format_mm9_refgene)
    print mm10
    mm10.save("hg38_refGene.glb")
    print mm10._findByLabel("name", "Nanog")