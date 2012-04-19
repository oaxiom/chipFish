#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display

"""

import sys, os, csv, time, cProfile, pstats

from .. import track
from .. import delayedlist # get delayedlist from glbase
from .. import genelist
from .. import format
from .. import config

def seqToTrk(infilename, outfilename, name, stranded=True, **kargs):
    """
    **Purpose**
        Convert a list of genomic coordinates to a trk database (Actually an SQL
        db)

    **Arguments**
        infilename
            the name of the filename to read from. This should be a csv.
            The data is loaded using a delayedlist(), so memory
            is not a big problem.

        outfilename
            The filename of the trk file to write to.

        name
            A name describing the track

        format (Optional, default=GERALD export.txt)
            a format specifier in the glbase format describing the csv
            file.

        stranded (Optional, default=True)
            expects a "strand" key and will store the strand.

    **Returns**
        True on completion
        and a trk file in outfilename

    """
    assert os.path.realpath(infilename), "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"

    if "format" in kargs and kargs["format"]:
        default_format = kargs["format"]
    else:
        default_format = format.exporttxt_loc_only

    if "force_cache" in kargs and kargs["force_cache"]: # This is undocumented and generally a really bad idea
        seqfile = genelist(filename=os.path.realpath(infilename), format=default_format)
    else:
        seqfile = delayedlist(filename=os.path.realpath(infilename), format=default_format)

    n = 0
    m = 0
    total = 0

    t = track(filename=outfilename, stranded=stranded, new=True, name=name, **kargs)

    config.log.info("Started %s -> %s" % (infilename, outfilename))

    s = time.time()
    for item in seqfile:
        if "strand" in item and stranded:
            t.add_location(item["loc"], strand=item["strand"])
        else:
            t.add_location(item["loc"])

        n += 1
        total += 1
        if n > 1000000:
            m += 1
            n = 0
            config.log.info("%s,000,000 tags read" % m) # How to correct for 1000 M tags?
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

    config.log.info("Finalise library. Contains '%s' tags" % total)
    t.finalise()
    config.log.info("Took: %s seconds" % (e-s))
    return(True)

if __name__ == "__main__":
    # testing:
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

