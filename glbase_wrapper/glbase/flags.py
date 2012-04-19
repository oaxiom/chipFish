"""
flags.py

. should be renamed helpers...
. This file is scheduled for deletion
"""

import utils, csv, cPickle

from helpers import strandSorter
from data import *

"""
valid accessory tags:

"any_tag": {"code": "code_insert_as_string"} # execute arbitrary code to construct this key.
"dialect": csv.excel_tab # dialect of the file, default = csv, set this to use tsv. or sniffer
"skip_lines": number # number of lines to skip at the head of the file.
"skiptill": skip until I see the first instance of <str>

"""

# lists of format-specifiers.
default = {"sniffer": 0} # the default, it loads your file based on the heading labels in the csv.
sniffer = default # alternative name
sniffer_tsv = {"sniffer": 0, "dialect": csv.excel_tab} # alternative name
