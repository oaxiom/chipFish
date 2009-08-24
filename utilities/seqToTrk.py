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


import sys, os, csv, cPickle

sys.path.append(os.path.realpath("../"))
from error import AssertionError
from track import track
from glbase_wrapper import delayedlist

def seqToTrk(infilenamepath, outfilenamepath):
    assert os.path.realpath(infilenamepath), "no filename specified"
    assert os.path.realpath(outfilenamepath), "no save filename specified"

    t = track(stranded=False, delayed=False) # strands not currently supported :(

    gerald_format = {"loc": {"code": "location(chr=column[1], left=column[12], right=int(column[12])+25)"}, "strand": 13,
        "dialect": csv.excel_tab}
    # strand is F/R ??

    seqfile = delayedlist(filename=os.path.realpath(infilenamepath), format=gerald_format)
    n = 0
    m = 0

    for item in seqfile:
        print item
        t.add_location(item["loc"])

        n += 1
        if n > 1000:
            m += 1
            n = 0
            print "%s,000,000" % m
            break

    oh = open(os.path.realpath(outfilenamepath), "wb")
    cPickle.dump(t, oh, -1)
    oh.close()

if __name__ == "__main__":
    # testing:
    print "Info: This may take a while..."
    seqToTrk("/home/hutchinsa/ChIP_Raw/CMN019_121_unique_hits.txt", "../data/out.trk")

    """
    print "Command-line seqToTrk:"
    print "Usage: $ ./seqToTrk.py <infile> <outfile>"

    assert sys.argv[1], "no infile specified"
    assert sys.argv[2], "no outfile specified"
    assert os.path.exists(sys.argv[1]), "infile does not exist"

    infile = sys.argv[1]
    outfile = sys.argv[2]
    """
