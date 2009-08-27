#! /usr/bin/python
"""

Utilties:
Part of chipFish,
converts a sequence file to a track (graph)-like display

Can also be used from the command-line, but only from this location
in chipFish. It's highly tied into chipFish and glbase
and not so easily stripped out.

(c) 2009 oAxiom

"""


import sys, os, csv, cPickle, time

sys.path.append(os.path.realpath("../"))
from error import AssertionError
from track import track
from glbase_wrapper import delayedlist

def seqToTrk(infilename, outfilename):
    assert os.path.realpath(infilename), "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"

    gerald_format = {"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13,
        "dialect": csv.excel_tab}
    # strand is F/R ??

    seqfile = delayedlist(filename=os.path.realpath(infilename), format=gerald_format)
    n = 0
    m = 0

    t = track(filename=outfilename, stranded=False, new=True) # strands not currently supported :(

    s = time.time()
    for item in seqfile:
        t.add_location(item["loc"], strand=item["strand"])

        n += 1
        if n > 1000000:
            m += 1
            n = 0
            print "%s,000,000" % m
    e = time.time()
    # 1000 = averages 8-9 s

    print "Took: %s s" % (e-s)

    t.finalise()
    t.show_debug_info()
    return(True)

if __name__ == "__main__":
    # testing:
    print "Info: This may take a while..."
    seqToTrk("/home/hutchinsa/ChIP_Raw/CMN019_121_unique_hits.txt", "../data/NSMash1.trk")

    """
    print "Command-line seqToTrk:"
    print "Usage: $ ./seqToTrk.py <infile> <outfile>"

    assert sys.argv[1], "no infile specified"
    assert sys.argv[2], "no outfile specified"
    assert os.path.exists(sys.argv[1]), "infile does not exist"

    infile = sys.argv[1]
    outfile = sys.argv[2]
    """
