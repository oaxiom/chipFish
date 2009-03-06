"""
utils, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for utilities for conversion etc...

"""

def rc(seq):
    """
    get the reverse complemnt of seq, returns rc_seq
    """
    compdict = {'A': 'T',
                'C': 'G',
                'G': 'C',
                'T': 'A',
                'N': 'N'
                }
    tseq = [] # new list
    for i in range(len(seq)):
        tseq.append(compdict[seq[i]])

    tseq.reverse()

    rc_seq = str(tseq).strip('[]').replace('\'', '').replace(' ', '').replace(',', '')
    return(rc_seq)

def getLocation(_string):
    """
    # split the string into a dict [chr, left, right]
    # string of for chrN:left-right
    """
    try:
        t = _string.split(":")
    except ValueError:
        return(False)

    chr = t[0].lower().strip("chr")
    if chr == "x":
        chr = "X"
    elif chr == "y":
        chr = "Y"
    elif chr == "m":
        chr = "M"
    else:
        chr = int(chr)

    t2 = t[1].split("-")
    left = int(t2[0])
    right = int(t2[1])
    return({"chr": chr, "left":left, "right":right})
