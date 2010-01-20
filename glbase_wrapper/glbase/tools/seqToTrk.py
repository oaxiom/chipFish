#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display


"""

import sys, os, csv, time, cProfile, pstats

sys.path.append(os.path.realpath("../"))
from error import AssertionError # these imports are why it can't be moved...
from track import track # the chipFish track implementation.
import delayedlist # get delayedlist from glbase

def seqToTrk(infilename, outfilename, **kargs):
    assert os.path.realpath(infilename), "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"

    gerald_format = {"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13,
        "dialect": csv.excel_tab}
    # strand is F/R ??
    #mikkelsen_format = {"loc": {"code": "location(chr=column[0].strip(\".fa\"), left=column[1], right=int(column[2]))"},
    #    "strand": 3,
    #    "dialect": csv.excel_tab} # this is actually a bed?

    seqfile = delayedlist(filename=os.path.realpath(infilename), format=gerald_format)
    n = 0
    m = 0

    t = track(filename=outfilename, stranded=False, new=True, **kargs) # strands not currently supported :(

    s = time.time()
    for item in seqfile:
        t.add_location(item["loc"], strand=item["strand"])

        n += 1
        if n > 1000000:
            m += 1
            n = 0
            print "%s,000,000" % m
            #break
    e = time.time()
    # 1000 = averages 8-9 s
    # 3.65 s cache.
    # 0.61 s better cacheing, less commits
    # 10,000 used to check now.
    # 5.98 s, 22.8Mb (Binary(array.tostring())
    # 6.0 s 17.1 Mb (str/eval)
    # 6.0 s 259kb (zlib/Binary()/tostring())
    # track_new
    # 2.5 s # overhead is almost all delayedlist now...

    print "Finalise:"
    t.finalise()
    print "Took: %s s" % (e-s)
    return(True)

if __name__ == "__main__":
    # testing:
    print "Info: This may take a while..."
    PROFILE = False
    if PROFILE:
        cProfile.run("seqToTrk(\"/home/hutchinsa/ChIP_Raw/CMN019_121_unique_hits.txt\", \"../data/NSMash1.trk\")", "seqToTrk.profile")

        p = pstats.Stats( "seqToTrk.profile")

        p.strip_dirs().sort_stats('time').print_stats()

    else:
        print "Command-line seqToTrk:"
        print "Usage: $ ./seqToTrk.py <infile> <outfile>"

        assert sys.argv[1], "no infile specified"
        assert sys.argv[2], "no outfile specified"
        assert os.path.exists(sys.argv[1]), "infile does not exist"

        infile = sys.argv[1]
        outfile = sys.argv[2]
