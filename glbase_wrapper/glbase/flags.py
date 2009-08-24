"""
flags.py

. should be renamed helpers...
"""

import utils, csv, cPickle

from helpers import strandSorter, barSplitter
from location import location
from data import ignorekeys, typical_headers

positive_strand_labels = frozenset(["+", "0", "f", "F"]) # correct? Sometimes 0 is the top strand?
negative_strand_labels = frozenset(["-", "1", "r", "R"])

# Template for doc strings.
"""
splurge

**Arguments**
    arg
        splurger

**Result**
    splurge
"""

regex_dict = {
    "a" : "a",
    "c" : "c",
    "g" : "g",
    "t" : "t",
    "r" : "[ag]", # twos
    "y" : "[ct]",
    "k" : "[gt]",
    "m" : "[ac]",
    "s" : "[gc]",
    "w" : "[at]",
    "h" : "[act]", # threes
    "v" : "[acg]",
    "d" : "[agt]",
    "b" : "[cgt]",
    "n" : "[acgt]" # four
}

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

# standard lists:
format_bed = {"loc": {"code": "utils.coordsToLocation(column[0], column[1], column[2])"}, "strand": 4, "dialect": csv.excel_tab,
    "skiplines": -1}

format_fasta = {"special": "fasta"}

# a selection of format files.
format_ann_list = {"loc": 4, "strand": 6, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1}# output format from my annotation script
format_peak_loc = {"loc": 4, "name": 9, "refseq": 11, "entrez": 12, "tag_height": 1} # A list of peak locations, annotated?
format_gs_out = {"refseq": 4, "entrez": 8, "fold_change": 1, "array_systematic_name": 0, "dialect": csv.excel_tab}# A list out of GeneSpring
format_array_simplified = {"refseq": 10, "array_systematic_name": 0, "entrez": 8}
format_illumina_anotations = {"array_systematic_name":0, "refseq": 3, "entrez": 7}

# load in MACS files:
format_macs_peak_loc = {"loc": 0, "tag_height": {"code": "int(column[4])", "fold": 8}} # format for coord modified macs file.

format_macs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 5,
    "skiplines": 13, "dialect": csv.excel_tab}

# Load in CCAT files:
format_ccat_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "dialect": csv.excel_tab,
    "tag_height": {"code": "int(barSplitter(column[3])[0])"}, "skiplines": -1}

format_ccat_output_csv = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"},
    "tag_height": {"code": "int(barSplitter(column[3])[0])"}, "fold_change": {"code": "barSplitter(3)[2]"}, "skiplines": -1}

# load in SISSRS file.
format_sissrs_output = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"}, "tag_height": 3,
    "fold": 4, "p-value": 5, "dialect": csv.excel_tab, "skiplines": 57}

# this is the format for array:
format_array_tsv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False,
        "dialect": csv.excel_tab}

format_array_csv = {"refseq": 0, "entrez": 1, "symbol": 2,
        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False}

# this is for genome: loads an intersect of knownGeneList vs refseq
format_genomeA = {"loc": "location(chr=column[0], column[2], column[3])",
        "strand": 1, "name": 2, "description": 3,
        "refseq": 4, "entrez": 5, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"}}

format_genomeB = {"loc": {"code": "location(chr=column[0], left=column[2], right=column[3])"}, "strand": 1, "name": 4,
    "description": 5, "refseq": 6, "entrez": 7,"dialect": csv.excel_tab,
    "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"}}

format_genomeC = {"loc": {"code": "location(chr=column[0], left=column[2], right=column[3])"},
        "strand": 1, "name": 2, "dialect": csv.excel_tab,
        "refseq": 3, "entrez": 4, "tss_loc": {"code": "strandSorter(column[0], column[2], column[3], column[1])"}}

# hg18 default refseq export.
format_hg18_refseq = {"loc": {"code": "location(chr=column[0], left=column[1], right=column[2])"},
        "strand": 5, "dialect": csv.excel_tab,
        "refseq": 3, "tss_loc": {"code": "strandSorter(column[0], column[1], column[2], column[5])"}}
