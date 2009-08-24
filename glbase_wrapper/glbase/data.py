"""

**data.py**

data stuff, part of glbase

Contains various data that doens't belong in helpers.py, flags.py or opt.py.
And probably shouldn't be exported outside of glbase
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

compdict = {'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'N': 'N',
            "a": "t",
            "c": "g",
            "g": "c",
            "t": "a",
            "n": "n"
            }

ignorekeys = frozenset(["dialect", "duplicates_key",
"skiplines", "debug", "special"]) # these are functional tags - so I should ignore them.

typical_headers = frozenset(["chipseq_loc", "loc", "chr", "#",
"Gene Name", "", "GenBank", "RefSeq",
"Systematic", "mm8.refGene.chrom", "mm8", "loc",
"mm9.refGene.chrom", "mm9",
"======================================================================", # stupid sissrs format garbage.
"=====================================================================", # stupid sissrs format garbage.
"======================================================================="] # stupid sissrs format garbage.
) # typical header labels;
