"""
utils, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for utilities for conversion etc...
stuff that doesn't really belong attached to a class.
useful macros that get used all over the code base

Later compile this one with Cython? Or split into a cutils.py for
macros needed fast?
"""

def rc(seq):
    """
    get the reverse complemnt of seq, returns rc_seq
    (simple, only deals with actgn)
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

def getLocation(location):
    """
    # split the string into a dict [chr, left, right]
    # string of for chrN:left-right
    returns string for chr and integers for left and right.
    """
    try:
        t = location.split(":")

        chr = t[0].lower().strip("chr")

        t2 = t[1].split("-")
        left = int(t2[0])
        right = int(t2[1])
        return({"chr": chr, "left":left, "right":right})
    except:
        return(False)

def formatLocation(chr, left, right):
    """
    returns a correctly formatted genomic location string.
    the opposite of getLocation()
    """
    return("chr%s:%s-%s" % (str(chr), str(left), str(right)))

def collide(leftA, rightA, leftB, rightB):
    """
    test collision between two 2D like objects.
    Candidate for C implementation.
    assumes left < right.
    """
    if leftA > rightB: return(False) # quickest failures.
    if rightA < leftB: return(False)

    if leftA > leftB and leftA < rightB:
        return(True)

    if rightA > leftB and rightA < rightB:
        return(True)

    if leftA < leftB and rightA > rightB:
        return(True)

    return(False)

if __name__ == "__main__":
    import time

    print collide(10, 20, 2, 8) # False, left of A.
    print collide(10, 20, 30, 40) # False, right of A.
    print collide(10, 20, 1, 30) # True, B extends over A
    print collide(10, 20, 15, 30) # True B within A
    print collide(10, 20, 5, 15) # True B within A.
    print collide(10, 20, 11, 18) # True B wholy within A.

    # takes 6 seconds to do 6 million collision checks...
    # too slow...

    s = time.time()

    for n in xrange(1000000):
        # all five conditions.
        a = collide(10, 20, 2, 8) # False, left of A.
        b = collide(10, 20, 30, 40) # False, right of A.
        c = collide(10, 20, 1, 30) # True, B extends over A
        d = collide(10, 20, 15, 30) # True B within A
        e = collide(10, 20, 5, 15) # True B within A.
        f = collide(10, 20, 11, 18) # True B wholy within A.

    e = time.time()
    print "time:", e - s
