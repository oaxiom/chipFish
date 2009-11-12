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


import sys, os, csv, time, cProfile, pstats

sys.path.append(os.path.realpath("../"))
from error import AssertionError # these imports are why it can't be moved...
from track import track # the chipFish track implementation.
from glbase_wrapper import delayedlist # get delayedlist from glbase

def seqToTrk(infilename, outfilename, **kargs):
    assert os.path.realpath(infilename), "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"

    gerald_format = {"loc": {"code": "location(chr=column[10].strip(\".fa\"), left=column[12], right=int(column[12])+25)"},
        "strand": 13,
        "dialect": csv.excel_tab}
    # strand is F/R ??
    mikkelsen_format = {"loc": {"code": "location(chr=column[0].strip(\".fa\"), left=column[1], right=int(column[2]))"},
        "strand": 3,
        "dialect": csv.excel_tab} # this is actually a bed?

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
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CMN019_121_unique_hits.txt", "../data/NSMash1_new.trk", name="NS5 cells Mash1")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME030_154_unique_hits.txt", "../data/SpMash1_new.trk", name="Spinal Cord Mash1")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/SCS843_053_reprocessed_unique_hits.txt", "../data/TcMash1_new.trk", name="Telencephalon Mash1")
        seqToTrk("/home/hutchinsa/ChIP_Raw/CMN010_010_reprocessed_unique_hits.txt", "../data/NSGFP_new.trk")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME031_031_unique_hits.txt", "../data/SpInput_new.trk")
        seqToTrk("/home/hutchinsa/ChIP_Raw/SCS844_844_reprocessed_unique_hits.txt", "../data/TcInput_new.trk")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME023_095_reprocessed_unique_hits.txt", "../data/ESRcor1.trk")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME024_096_reprocessed_unique_hits.txt", "../data/ESRcor2.trk")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME025_097_reprocessed_unique_hits.txt", "../data/ESRcor3.trk")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME026_099_reprocessed_unique_hits.txt", "../data/ESRest.trk")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME027_100_reprocessed_unique_hits.txt", "../data/ESSin3a.trk")
        #seqToTrk("/home/hutchinsa/ChIP_Raw/CME028_101_reprocessed_unique_hits.txt", "../data/ESSin3b.trk")
        # The Mikkelsen histone ChIP-seq data:
        # the name attribute is not actualyl implemented...
        """
        seqToTrk("/home/hutchinsa/Histone marks/NP_K4.txt", "../data/NS_H3K4me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/NP_K27.txt", "../data/NS_H3K27me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/NP_K36.txt", "../data/NS_H3K36me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/NP_K9.txt", "../data/NS_H3K9me3.trk", name="NS5 GFP Control")

        seqToTrk("/home/hutchinsa/Histone marks/ES_K4.txt", "../data/ES_H3K4me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/ES_K27.txt", "../data/ES_H3K27me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/ES_K36.txt", "../data/ES_H3K36me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/ES_K9.txt", "../data/ES_H3K9me3.trk", name="NS5 GFP Control")

        seqToTrk("/home/hutchinsa/Histone marks/MEF_K4.txt", "../data/MEF_H3K4me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/MEF_K27.txt", "../data/MEF_H3K27me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/MEF_K36.txt", "../data/MEF_H3K36me3.trk", name="NS5 GFP Control")
        seqToTrk("/home/hutchinsa/Histone marks/MEF_K9.txt", "../data/MEF_H3K9me3.trk", name="NS5 GFP Control")
        """

    """
    print "Command-line seqToTrk:"
    print "Usage: $ ./seqToTrk.py <infile> <outfile>"

    assert sys.argv[1], "no infile specified"
    assert sys.argv[2], "no outfile specified"
    assert os.path.exists(sys.argv[1]), "infile does not exist"

    infile = sys.argv[1]
    outfile = sys.argv[2]
    """
