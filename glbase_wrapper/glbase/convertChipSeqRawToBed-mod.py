"""
"""

import sys, os, csv

def collide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
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

def convertRawToBed(path, _in, _out):
    print "convertRawToBed(path=%s _in=%s, _out=%s)" % (path, _in, _out)
    ih = open(os.path.join(path, _in), "rb")
    oh = open(os.path.join(path, _out), "wb")

    reader = csv.reader(ih, dialect=csv.excel_tab)
    writer = csv.writer(oh, dialect=csv.excel_tab)

    for line in reader:
        t = line[10].split(".")
        chr = t[0]
        left = int(line[12])
        if line[13] == "F":
            strand = "+"
        elif line[13] == "R":
            strand = "-"
        out = [chr, left, left+35, "Read", line[15], strand]
        writer.writerow(out)

    # out format, tsv:
    # [chr, left, left+36, "Read", ?, srand (+/-)
    ih.close()
    oh.close()

def getLocation(_string):
    # split the string into a dict [chr, left, right]
    # string of for chrN:left-right
    # incorrectly returns the X and Y tags?!
    try:
        t = _string.split(":")
    except ValueError:
        return(False)

    chr = t[0].lower().strip("chr")
    t2 = t[1].split("-")
    left = int(t2[0])
    right = int(t2[1])
    return({"chr": chr, "left":left, "right":right})

def findChipSeqHit(path, _in, _out, coords):
    """
    Scan an export.txt file and output all sequencing hits within the coords.
    puts it into a format suitable to upload to UCSC to view.
    """
    print "findChipSeqHit(path=%s, _in=%s, _out=%s, coords=%s)" % (path, _in, _out, coords)
    ih = open(os.path.join(path, _in), "rb")
    oh = open(os.path.join(path, _out), "wb")

    reader = csv.reader(ih, dialect=csv.excel_tab)
    writer = csv.writer(oh, dialect=csv.excel_tab)

    n = 0
    m = 0
    f = 0
    range = getLocation(coords)
    for line in reader:
        t = line[10].split(".")
        chr = t[0].strip("chr")
        loc = int(line[12])
        n += 1
        if n >1000000:
            m += 1
            print m,"million - found:", f
            n = 0

        if chr == range["chr"]:
            if collide(range["left"], range["right"], loc, loc+35):
                chr = line[10].split(".")[0]
                left = int(line[12])
                if line[13] == "F":
                    strand = "+"
                    col = "255,0,0"
                    writer.writerow([chr, left, left+135, "read", line[15], strand, left, left+35, col])
                elif line[13] == "R":
                    strand = "-"
                    col = "0,0,255"
                    writer.writerow([chr, left-100, left+35, "read", line[15], strand, left, left+35, col])
                f += 1

    print "Found:", f

    oh.close()
    ih.close()

def topXlines(path, file, X):
    # print the top X lines of a csv file;
    oh = open(os.path.join(path, file), "rb")
    reader = csv.reader(oh, dialect=csv.excel_tab)

    read = 0
    for line in reader:
        print line

        read += 1
        if read > X:
            break
    oh.close()

def getAllSequencesInRange(path, _in, _out, location):
    """
    Scan an export.txt file and output all sequencing hits within the coords.
    """
    print "getAllSequencesInRange(path=",path," _in=",_in," _out=", _out, location,")"
    ih = open(os.path.join(path, _in), "rb")
    oh = open(os.path.join(path, _out), "wb")

    reader = csv.reader(ih, dialect=csv.excel_tab)
    writer = csv.writer(oh)

    n = 0
    m = 0
    range = getLocation(location)
    for line in reader:
        t = line[10].split(".")
        chr = t[0].strip("chr")
        loc = int(line[12])
        n += 1
        if n >1000000:
            m += 1
            print m,"million"
            n = 0

        if str(chr) == str(range["chr"]):
            if collide(range["left"], range["right"], loc, loc+35):
                writer.writerow(line)
                print "Found: 1"

    oh.close()
    ih.close()

def convertRawToALN(path, _in, _out):
    """
    For converting into a format suitable for cisgenome (it's native ALN file?)
    """
    print "convertRawToALN(path=",path," _in=",_in," _out=", _out,")"
    ih = open(os.path.join(path, _in), "rb")
    oh = open(os.path.join(path, _out), "wb")

    reader = csv.reader(ih, dialect=csv.excel_tab)
    writer = csv.writer(oh, dialect=csv.excel_tab)

    for line in reader:
        t = line[10].split(".")
        chr = t[0]
        left = int(line[12])
        out = [chr, left, line[13]]
        writer.writerow(out)

    ih.close()
    oh.close()

if __name__ == "__main__":
    path = "/home/hutchinsa/ChIP_Raw/"
    #convertRawToBed(path, "Top_10_lines.txt", "Top_10_lines_mine.bed")
    #convertRawToBed(path, "SCS845_062_reprocessed_unique_hits.txt", "SCS845_re_unique_hits.bed") # NSSox2
    #convertRawToBed(path, "SCS843_053_reprocessed_unique_hits.txt", "SCS843_re_unique_hits.bed")
    #convertRawToBed(path, "SCS846_846_reprocessed_unique_hits.txt", "SCS846_re_unique_hits.bed") # F9GFP
    #convertRawToBed(path, "SCS847_108_reprocessed_unique_hits.txt", "SCS847_re_unique_hits.bed") # F9Oct4
    #convertRawToBed(path, "SCS848_109_reprocessed_unique_hits.txt", "SCS848_re_unique_hits.bed") # F9Sox17
    #convertRawToBed(path, "CMN010_010_reprocessed_unique_hits.txt", "CMN010_re_unique_hits.bed") # NSControl (not ours though)
    #convertRawToBed(path, "CMN011_084_reprocessed_unique_hits.txt", "CMN011_re_unique_hits.bed") # NSMash data
    convertRawToBed(path, "CMN019_121_unique_hits.txt", "CMN019_121_unique_hits.bed") # NSMash data
    #topXlines(path, "CMN011_084_reprocessed_unique_hits.txt", 5)
    topXlines(path, "/home/hutchinsa/Histone marks/NP_WCE.txt", 5)
    topXlines(path, "SCS845_re_unique_hits.bed", 5)

    #topXlines(path, "SCS846_846_unique_hits.bed", 1)
    #topXlines(path, "Top_10_lines.txt", 1)
    #convertRawToALN(path, "SCS845_062_unique_hits.txt", "SCS845_unique_hits.aln")
    #topXlines(path, "SCS845_unique_hits.aln", 10)
    #convertRawToALN(path, "CMN011_084_unique_hits.txt", "CMN011_unique_hits.aln")
    #getAllSequencesInRange(path, "SCS845_062_reprocessed_unique_hits.txt", "NSSox2_Nestin_sites.csv", "chr3:88060027-88060776")
    #findChipSeqHit(path, "SCS845_062_reprocessed_unique_hits.txt", "NSSox2_Nestin_sites.bed", "chr3:88060027-88060776")
    #findChipSeqHit(path, "CMN010_010_reprocessed_unique_hits.txt", "NSCon_Nestin_sites.bed", "chr3:88060027-88060776")
    #findChipSeqHit(path, "SCS845_062_reprocessed_unique_hits.txt", "NSSox2_Olig1_sites.bed", "chr16:91151613-91152131")
    #findChipSeqHit(path, "CMN010_010_reprocessed_unique_hits.txt", "NSCon_Olig1_sites.bed", "chr16:91151613-91152131")
    #findChipSeqHit(path, "SCS845_062_reprocessed_unique_hits.txt", "NSSox2_Sox2_sites.bed", "chr3:34845391-34845817")
    #findChipSeqHit(path, "CMN010_010_reprocessed_unique_hits.txt", "NSCon_Sox2_sites.bed", "chr3:34845391-34845817")
    #getAllSequencesInRange(path, "Top_10_lines.txt", "NSSox2_Nestin_sites.csv", "chr10:88060127-88060676")
    #findChipSeqHit(path, "CMN011_084_reprocessed_unique_hits.txt", "NSMash1_Nestin_sites.bed", "chr3:88060027-88060776")
    #findChipSeqHit(path, "SCS845_062_reprocessed_unique_hits.txt", "NSSox2_Olig1_sites.bed", "chr16:91151613-91152131")
    #findChipSeqHit(path, "CMN011_084_reprocessed_unique_hits.txt", "NSMash1_Olig1_sites.bed", "chr16:91151613-91152131")
    #findChipSeqHit(path, "SCS845_062_reprocessed_unique_hits.txt", "NSSox2_Sox2_sites.bed", "chr3:34845391-34845817")
    #findChipSeqHit(path, "CMN011_084_reprocessed_unique_hits.txt", "NSMash1_Sox2_sites.bed", "chr3:34845391-34845817")

    # EmMash1 data
    #findChipSeqHit(path, "SCS843_053_reprocessed_unique_hits.txt", "EmMash1_lfng_sites.bed", "chr5:140863497-140863902")
    #findChipSeqHit(path, "SCS843_053_reprocessed_unique_hits.txt", "EmMash1_Dll1_A_sites.bed", "chr17:15081928-15082138")
    #findChipSeqHit(path, "SCS843_053_reprocessed_unique_hits.txt", "EmMash1_Dll1_B_sites.bed", "chr17:15081928-15082138")
    #chr17:15081897-15082143
    #chr17:15076967-15077091
    """
    findChipSeqHit(path, "SCS848_109_reprocessed_unique_hits.txt", "F9Sox17_LmnA_1_sites.bed", "chr3:88579026-88580226")
    findChipSeqHit(path, "SCS848_109_reprocessed_unique_hits.txt", "F9Sox17_LmnA_2_sites.bed", "chr3:88589059-88590259")
    findChipSeqHit(path, "SCS848_109_reprocessed_unique_hits.txt", "F9Sox17_LmnA_3_sites.bed", "chr3:88589059-88590259")

    findChipSeqHit(path, "SCS847_108_reprocessed_unique_hits.txt", "F9Oct4_LmnA_1_sites.bed", "chr3:88579026-88580226")
    findChipSeqHit(path, "SCS847_108_reprocessed_unique_hits.txt", "F9Oct4_LmnA_2_sites.bed", "chr3:88589059-88590259")
    findChipSeqHit(path, "SCS847_108_reprocessed_unique_hits.txt", "F9Oct4_LmnA_3_sites.bed", "chr3:88589059-88590259")
    """
