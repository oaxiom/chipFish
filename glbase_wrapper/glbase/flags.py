"""
flags.py

. should be renamed helpers...
"""

import utils, csv, cPickle

from helpers import strandSorter
from data import *

"""
valid accessory tags:

"any_tag": {"code": "code_insert_as_string"} # execute arbitrary code to construct this key.
"dialect": csv.excel_tab # dialect of the file, default = csv, set this to use tsv. or sniffer
"skip_lines": number # number of lines to skip at the head of the file.

not yet implemneted, but suggested:
"wait_for_line": "a string to wait for in column[0]"

"""

# lists of format-specifiers.
default = {"sniffer": 0} # the default, it loads your file based on the heading labels in the csv.
sniffer = default # alternative name
sniffer_tsv = {"sniffer": 0, "dialect": csv.excel_tab} # alternative name

# standard lists:
format_bed = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "strand": 4, "dialect": csv.excel_tab,
    "skiplines": -1}
format_minimal_bed = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "dialect": csv.excel_tab,
    "skiplines": -1} # no strand basically, and ignore any other keys.
format_bed_no_strand = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "dialect": csv.excel_tab,
    "skiplines": -1}

format_fasta = {"special": "fasta"}

# a selection of format files.
format_ann_list = {"loc": 4, "strand": 6, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1}# output format from my annotation script
format_peak_loc = {"loc": 4, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1} # A list of peak locations, annotated?
format_gs_out = {"refseq": 4, "entrez": 8, "fold_change": 1, "array_systematic_name": 0, "dialect": csv.excel_tab}# A list out of GeneSpring
format_array_simplified = {"refseq": 10, "array_systematic_name": 0, "entrez": 8}
format_illumina_anotations = {"array_systematic_name":0, "refseq": 3, "entrez": 7}

format_gtf = {"feature_type": 1, "feature": 2, "gtf_decorators": 8,
	"loc": {"code": "location(chr=column[0], left=column[3], right=column[4])"},  
	"skiplines": -1, "dialect": csv.excel_tab, "debug": True}

# load in MACS files:
format_macs_peak_loc = {"loc": 0, "tag_height": {"code": "int(column[4])", "fold": 8}} # format for coord modified macs file.

format_macs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 5,
    "skiptill": "chr", "dialect": csv.excel_tab, "fold_change": 7}

exporttxt_loc_only_format = {"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13,
        "dialect": csv.excel_tab} # export.txt file (output from the illumina pipeline), but only loads the location and strand.

exporttxt_all_format = {"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13, "seq": 6, "quality_score": 7,
        "dialect": csv.excel_tab} # export.txt file (output from the illumina pipeline), but only loads the location and strand.

# This is not a full implementation of the sam file specification.
# but it will accept a sam file as outputted by tophat.
sam_tophat = {"name": 0, "loc": {"code": "location(chr=column[2], left=column[3], right=column[3])"},
    "mapq": 3, "seq": 9, "dialect": csv.excel_tab, "skiplines": -1}

# Load in CCAT files:
# CCAT1.3 file format:
# <chromosome> <position of the peak> <start of region> <end of region> <read counts in ChIP library> <read counts in re-sampled control library> <the fold-change against control>

format_ccat_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[1])"}, "dialect": csv.excel_tab,
    "tag_height": 4, "fold": 6, "skiplines": -1}

format_ccat_output_csv = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[1])"},
    "tag_height": 4, "fold": 6, "fold_change": {"code": "barSplitter(3)[2]"}, "skiplines": -1}

# load in SISSRS file.
format_sissrs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 3,
    "fold": 4, "p-value": 5, "dialect": csv.excel_tab, "skiplines": 57}

# this is the format for array:
format_array_tsv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False,
        "dialect": csv.excel_tab}

format_array_csv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False}

# mm8 refGene table from UCSC:
format_mm8_refgene = {"loc": {"code": "location(chr=column[0], left=column[2], right=column[3])"},
        "strand": 1, "name": 4, "description": 5, "dialect" : csv.excel_tab,
        "refseq": 6, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"}
        }

# mm9 refGene table from UCSC:
format_mm9_refgene = {"loc": {"code": "location(chr=column[2], left=column[4], right=column[5])"},
        "strand": 3, "name": 12, "dialect" : csv.excel_tab,
        "refseq": 1, "tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"},
        "cds_loc": {"code": "location(chr=column[2], left=column[6], right=column[7])"},
        "exons_count": 8,
        } # description is lost from the mm9 table?

format_ensembl = {"loc": {"code": "location(chr=column[2], left=column[4], right=column[5])"},
	"tss_loc": {"code": "strandSorter(column[2], column[4], column[5], column[3])"},
	"ensmbl": 1, "name": 12, "exon_count": 8,
	"dialect": csv.excel_tab, "skiplines": 0} 
	
# hg18 default refseq export.
format_hg18_refseq = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"},
        "strand": 5, "dialect": csv.excel_tab,
        "refseq": 3, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"}}
