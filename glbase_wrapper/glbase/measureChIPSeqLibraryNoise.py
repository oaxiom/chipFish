"""
measureChIPSeqLibraryNoise.py
"""

import sys, os, csv, genelist, utils

def collide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    I really should turn this into a C library.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return(False)
    if Aleft > Bright:
        return(False)

    if Aleft <= Bright and Aright >= Bright:
        return(True) # Bright point is within A, thus collision

    if Aright >= Bleft and Aleft <= Bleft:
        return(True) # Bleft point is within A, thus collision.

    if Bleft <= Aright and Bright >= Aright:
        return(True) # Aright point is within B, thus collision

    if Bright >= Aleft and Bleft <= Aleft:
        return(True) # Aleft point is within B, thus collision.

    #print "unhandled!"
    return(False)

def do(path, chip_seq_filename, geneList):
    oh = open(os.path.join(path,  chip_seq_filename), "rU")
    reader = csv.reader(oh, dialect=csv.excel_tab)

    total = 0
    found = 0
    n = 1

    oD = geneList.getOrderedData()

    for item in reader:
        chr = item[0].strip("chr")
        if oD.has_key(chr):
            left = int(item[1])
            right = int(item[2])
            for tag in oD[chr]:
                tag_loc = utils.getLocation(tag["data"]["loc"])
                mid_tag = tag_loc["left"] - ((tag_loc["left"] - tag_loc["right"]) / 2)
                #print tag_loc, mid_tag
                if collide(mid_tag-200, mid_tag+200, left, right):
                    found += 1
                    #break
        total += 1
        n += 1
        if n > 1000000: # counter
            print "Done: %s - Found: %s tags, %f, " % (total, found, float(found)/total)
            n = 1
            #break

    print "Found:%s Total:%s Percent:%s Peaks:%s Tags/Peak:%s nomalised score:%s" % (found, total, (float(found) / total)*100.0, len(geneList), float(found)/len(geneList), ((float(found)*100.0/len(geneList)*100.0) / total))

if __name__ == "__main__":
    format_macs_peak_loc = {"loc": 0, "tag_height": 6, "fold": 8} # format for coord modified macs file.

    print "go!"

    oh = open("belh.txt", "wb")
    sys.stdout = oh

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "F9Oct4_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "SCS847_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "F9Sox17_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "SCS848_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "EmMash1_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "SCS843_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "NSSox2_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "SCS845_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "NSSox2_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "CMN011_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "NSMash1_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name, "correct"
    do(path, "CMN011_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "NSMash1_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "CMN019_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "NSMash1_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "CMN010_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "MSMash1_WCL_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name, "correct"
    do(path, "CMN019_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "MSMash1_WCL_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "CMN011_re_unique_hits.bed", g)

    path = "/home/hutchinsa/ChIP_Raw/"
    g = genelist.geneList(path, "MSMash1_WCL_peaks.csv", format_macs_peak_loc, "loc", True)
    print g.name
    do(path, "CMN010_re_unique_hits.bed", g)
