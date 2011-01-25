"""
Utilities

Various utilities to support the genome scanning scripts

R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
degenerate character N=[ATCG]
3-fold degenerate motifs re not used like the Lander paper.

"""

import sys, os, numpy, string, csv, random, math

#from errors import AssertionError

def library(args):
    """
    lovely generator...
    I'm still not sure exactly how this works...
    """
    if not args:
        yield ""
        return
    for i in args[0]:
        for tmp in library(args[1:]):
            yield i + tmp
    return

def expandDegenerateMotifs(motif):
    """
    expand a degeneratemotif into its constituant parts;
    returns a list;
    motif should be of the form:
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    degenerate character N=[ATCG]
    e.g. "raangt"
    """
    mlen = len(motif)
    nm = motif.lower() # just in case;

    # scan the motif and count the no of degenerate motifs
    newlist = []
    for n in xrange(mlen):
        if nm[n] == "r" or nm[n] == "y" or nm[n] == "k" or nm[n] == "m" or nm[n] == "s" or nm[n] == "w": # AG
            newlist.append((2, n, nm[n])) # append a triple, with the number of flips and its location
        elif nm[n] == "n":
            newlist.append((4, n, nm[n]))
        else:
            newlist.append((1, n, nm[n]))

    if len(newlist) == 0 : return([motif]) # no degenerate elements in the motif

    l = []
    #print newlist
    iti(newlist, 0, None, l)

    return(l)

def iti(_fm, _fmcpos, _cl, _list): # my iterator
    """
    This is possibly the best piece of code I've ever made!
    Mainly because it's almost entirely unitelligable....
    It's some sort of iterator to to generate
    N-mers.
    I think this has been replaced by regexs now.
    fm = a triple of the form (number to flip, location, seq_at_pos)
    """
    # do some set up
    if not _cl: # also okay if _cl == False; be careful with these, as may be False, but not None
        _cl = []
        for n in xrange(len(_fm)):
            _cl.append("")

    # the main iterator;
    if _fmcpos >= len(_fm):
        return(False) # reached end of list, signal to stick it back on the end;
    else:
        n = _fm[_fmcpos]
        #print n
        for x in xrange(n[0]):
            _cl[n[1]] = osc(_cl[n[1]], n[2])
            if not iti(_fm, _fmcpos+1, _cl, _list): # each time we iterate at the end of the motif add it to the list;
                # convert the list back to a string
                _copy = string.join(_cl, '')
                _list.append(_copy)
        return(True)
    return(True)

def movingAverage(listIn, window=20, normalise=False, bAbsiscaCorrect=True):
    assert window < len(listIn), "the window size for the moving average is too large"
    assert window >= 1, "the window size is too small (%s < 1)" % window

    if window == 1: # just return the original array
        return(numpy.arange(0, len(listIn)), listIn)


    if bAbsiscaCorrect:
        half_window_left = int(math.ceil(window / 2.0)) # correct for floating error division.
        half_window_right = int(math.floor(window / 2.0))
        x = numpy.arange(half_window_left, len(listIn)-half_window_right)
    else:
        x = numpy.arange(0, len(listIn)-window)

    y = []

    for n in xrange(half_window_left, len(listIn)-half_window_right):
        score = 0
        for i in xrange(n-half_window_left, n+half_window_right, 1):
            score += listIn[i]

        if normalise:
            y.append(float(score) / window)
        else:
            y.append(score)

    return(x, y)

def osc(last, type):
    """
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    degenerate character N=[ATCG]
    """
    if type == "r":
        if last == "a": return("g")
        if last == "g": return("a")
        return("a")
    if type == "y":
        if last == "c": return("t")
        if last == "t": return("c")
        return("c")
    if type == "k":
        if last == "g": return("t")
        if last == "t": return("g")
        return("g")
    if type == "m":
        if last == "a": return("c")
        if last == "c": return("a")
        return("a")
    if type == "s":
        if last == "g": return("c")
        if last == "c": return("g")
        return("g")
    if type == "w":
        if last == "a": return("t")
        if last == "t": return("a")
        return("a")
    if type == "n":
        if last == "a": return("c")
        if last == "c": return("g")
        if last == "g": return("t")
        if last == "t": return("a")
        return("a")
    return(type)

def rc(seq):
    """
    get the reverse complemnt of seq
    """
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
    tseq = [compdict[i] for i in seq] # new list
    tseq.reverse()

    return("".join(tseq))

def rc_expanded(seq):
    """
    get the reverse complemnt of seq
    this one works for the expanded alphabet;
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
    [ACT] = h, [ACG] = v, [AGT] = d, [CGT] = b
    degenerate character N=[ATCG]
    """
    compdict = {'a': 't',
                'c': 'g',
                'g': 'c',
                't': 'a',
                'r': 'y',
                'y': 'r',
                'k': 'm',
                'm': 'k',
                's': 's', # palindromic
                'w': 'w', # palindromic
                'n': 'n', # all;
                "h": "d", # threes
                "v": "b",
                "d": "h",
                "b": "v"
                }
    tseq = [compdict[i] for i in seq] # new list
    tseq.reverse()

    return("".join(tseq))

def expandElement(_element, _preserve=True): #
    """
    _element = (string) to expand

    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for left in xrange(4): # base iterator;
        for right in xrange(4):
            lib.append(dir[left]+_element+dir[right])

    return  lib

def expandElementRightOnly(_element, _preserve=True): #
    """
    _element = (string) to expand

    THis one only expands the element one base pair 3'

    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for right in xrange(4):
        lib.append(_element+dir[right])

    return  lib

def expandElementRightOnly_degenerate(_element, _preserve=True): #
    """
    _element = (string) to expand

    THis one only expands the element one base pair 3'
    will use a degenerate code;
    R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT], and the four-fold
degenerate character N=[ATCG]

    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t",
           4 : "r",
           5 : "y",
           6 : "k",
           7: "m",
           8: "s",
           9: "w",
           10:"n"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for right in xrange(11):
        lib.append(_element+dir[right])

    return  lib

def expandElementRightOnly_degenerate_n(_element, _preserve=True): #
    """
    _element = (string) to expand

    THis one only expands the element one base pair 3'
    will use a degenerate code;
    this one n only
    returns
    (list) of expanded elements
    """
    dir = {0 : "a",
           1 : "c",
           2 : "g",
           3 : "t",
           4:"n"}

    if _preserve:
        lib = [_element] # return string, preserve the previous element in lib[0]
    else:
        lib = []

    for right in xrange(5):
        lib.append(_element+dir[right])

    return  lib

def convertCSVtoFASTA(csvfilename, outfile, sequenceCol, fastaNameCol = None):
    """
    csvfilename = the csv filename to load without the path, assumed to be in environment.mouseGenomePath
    outfile = outfile name, also in environment.mouseGenomePath
    sequenceCol (integer) = the number of the column (starts from 0)
    fastaNameCol (integer) (optional) = the column to use as a FASTA name, otherwise will default to 0..n
    ignores a row if the first column begins with "#"
    """
    ofh = open(os.path.join(sys.path[0], csvfilename), "rb")
    sfh = open(os.path.join(sys.path[0], outfile), "wb")

    csvreader = csv.reader(ofh)

    fastanameseries = 0

    for n in csvreader:
        t = n[0]
        if t[0] <> "#":
            seq = n[sequenceCol]
            if fastaNameCol:
                fastname = n[fastaNameCol]
                sfh.write('>'+fastname+'\r\n')
            else: # no fasta name so use the series
                sfh.write('>'+str(fastanameseries)+'\r\n')
                fastanameseries += 1
            sfh.write(seq+'\r\n')

    ofh.close()
    sfh.close()

def convertFASTAtoCSV(filename):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    """
    assert os.path.exists(filename), "filename %s not found" % filename

    try:
        openfile = open(filename, "rb")
        savefile = open(filename+'_out.csv', "wb")
    except IOError:
        print "Error opening File"
        sys.exit()

    writer = csv.writer(savefile)

    record = ""
    entry = Node("empty")
    for line in openfile:
        if line[:1] != ">": # not a FASTA block, so add this line to the sequence
            entry.seq += line.replace('\r\n', '') # strip out the new line WINDOWS specific!

        if line[:1] == ">": # fasta start block
            # start recording
            # add the old Node to the list
            if entry.name != "empty":
                # convert the list to a tuple
                writer.writerow([entry.name, "", "", "", entry.seq])
                del entry
            entry = Node(line) # make a new node

def convertFASTAtoDict(filename):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    returns a list of the form [{name, seq}, ... {name, seq}]
    """
    assert os.path.exists(filename), "filename %s not found" % filename

    openfile = open(filename, "rU")

    result = []
    record = ""
    #entry = {"seq": "", "name": "empty"}
    for line in openfile:
        if line[:1] != ">": # not a FASTA block, so add this line to the sequence
            entry["seq"] += line.replace('\r', '').replace("\n", "") # strip out the new line WINDOWS specific!

        if line[:1] == ">": # fasta start block
            entry = {"seq": "", "name": "empty"} # make a new node
            # start recording
            entry["name"] = line.replace(">", "")
            # add the old Node to the list
            if entry["name"] != "empty":
                # convert the list to a tuple
                result.append(entry)
    return(result)

def scanNumberOfBasePairs(fastafilehandle):
    """
    pass me a fasta file handle (already open()'d) and this will return the raw count
    or some other iterable object;
    """
    dict = {"a" : 0,
            "A" : 0,
            "c" : 1,
            "C" : 1,
            "g" : 2,
            "G" : 2,
            "t" : 3,
            "T" : 3}

    a = []
    a.append(0)
    a.append(0)
    a.append(0)
    a.append(0)

    for line in fastafilehandle:
        if line[0] != ">": # there is a more elegant way to do this...
            lcline = line.lower()
            a[dict["a"]] += lcline.count("a")
            a[dict["c"]] += lcline.count("c")
            a[dict["g"]] += lcline.count("g")
            a[dict["t"]] += lcline.count("t")
    return(a)

def removeDuplicatesFromCSV(path, csvfile, outfile, column_no = 3, bKeepEmptyCols=False):
    """
    delete duplicates based on the column no
    """
    inf = open(os.path.join(path, csvfile), "rb")
    outf = open(os.path.join(path, outfile), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)
    ulist = []

    for line in reader:
        if line[column_no] in ulist:
            # don't write this enty,
            print "Duplicate:", line[column_no]
        else:
            # add to ulist and write to file;
            if line[column_no]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[column_no])
                writer.writerow(line) # if I tab this in - don't keep empty rows
            else: # col is empty;
                #print "Duplicate: <empty>"
                if bKeepEmptyCols: writer.writerow(line)
    inf.close()
    outf.close()

def collide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return(False)
    if Aleft > Bright:
        return(False)

    if Aleft == Bleft: return(1) # I have to cheat here otherwise it returns 0 which will evaluate as False;
    if Aleft == Bright: return(1)
    if Aright == Bleft: return(1)
    if Aright == Bright: return(1)

    if Aleft <= Bright and Aright >= Bright:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Bright point is within A, thus collision

    if Aright >= Bleft and Aleft <= Bleft:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Bleft point is within A, thus collision.

    if Bleft <= Aright and Bright >= Aright:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Aright point is within B, thus collision

    if Bright >= Aleft and Bleft <= Aleft:
        A = abs(Aleft - Bright)
        B = abs(Aright - Bright)
        C = abs(Aleft - Bleft)
        D = abs(Aright - Bleft)
        closest = min(A, B, C, D)
        return(closest) # Aleft point is within B, thus collision.

    #print "unhandled!"
    return(False)

def qcollide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return(False)
    if Aleft > Bright:
        return(False)

    if Aleft <= Bright and Aright >= Bright:
        return(True) # Bright point is within A, collision

    if Aright >= Bleft and Aleft <= Bleft:
        return(True) # Bleft point is within A, collision.

    if Bleft <= Aright and Bright >= Aright:
        return(True) # Aright point is within B, collision

    if Bright >= Aleft and Bleft <= Aleft:
        return(True) # Aleft point is within B, collision.

    #print "unhandled!"
    return(False)

def removeDuplicatesFromListOfDicts(list_of_dicts, key):
    """
    remove duplicates from a list of dicts based on key,
    returns the list in the format it arrived.
    only checks key. Does not check anything else.
    """
    ulist = []
    newlist = []
    dupecount = 0

    for line in list_of_dicts:
        if line[key] in ulist:
            # don't write this enty,
            dupecount +=1
        else:
            # add to ulist and write to file;
            if line[key]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[key])
            newlist.append(line)
    #print ">> Duplicates Found:", dupecount
    return(newlist)

def removeDuplicatesFrom2DList(_list, column_no = 3):
    """
    delete duplicates based on the column no
    returns a list
    deletes dupes in a 2D-csv like list;
    """
    ulist = []
    newlist = []
    dupecount = 0

    for line in _list:
        if line[column_no] in ulist:
            # don't write this enty,
            dupecount +=1
            print "Duplicate:%s" % (line[column_no])
        else:
            # add to ulist and write to file;
            if line[column_no]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[column_no])
            newlist.append(line)
    #print ">> Duplicates Found:", dupecount
    return(newlist)

def removeDuplicatesFromCSV_2Cols(path, csvfile, outfile, column_no1 = 0, column_no2 = 1):
    """
    delete duplicates based on the column no1 and 2
    """
    inf = open(os.path.join(path, csvfile), "rb")
    outf = open(os.path.join(path, outfile), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)
    ulist = []

    for line in reader:
        if line[column_no1]+line[column_no2] in ulist:
            # don't write this enty,
            print "duplicate:", line[column_no1]+line[column_no2]
        else:
            # add to ulist and write to file;
            ulist.append(line[column_no1]+line[column_no2])
            writer.writerow(line)
    inf.close()
    outf.close()

def keepRowOnlyIfColXHasValue(path, _in, _out, _colX):
    """
    what it says
    """
    print "keepRowOnlyIfColXHasValue(",path, _in, _out, _colX,")"
    inf = open(os.path.join(path, _in), "rb")
    outf = open(os.path.join(path, _out), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)

    for line in reader:
        if len(line) > _colX:
            if line[_colX]:
                writer.writerow(line)
    inf.close()
    outf.close()

def FASTAToCSV(filename):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    """
    #try:
    openfile = open(filename, "rb")
    savefile = open(filename+'_out.csv', "wb")
    #except IOError:
    #    print "Error opening File: %s" % filename
    #    sys.exit()

    writer = csv.writer(savefile)

    record = ""
    entry = Node("empty")
    for line in openfile:
        if line[:1] != ">": # not a FASTA block, so add this line to the sequence
            entry.seq += line.replace("\r", "").replace("\n", "")

        if line[:1] == ">": # fasta start block
            # start recording
            # add the old Node to the list
            if entry.name != "empty":
                # convert the list to a tuple
                writer.writerow([entry.name, "", "", entry.seq])
                del entry
            entry = Node(line) # make a new node

def FASTAToLIST(filename):
    """
    load a fasta file and output it into a big list;
    expects filename to be correct
    """
    try:
        openfile = open(filename, "rb")
    except IOError:
        print "Error opening File:", filename
        sys.exit()

    record = ""
    elementList = []
    entry = Node("empty")
    for line in openfile:
        if line[:1] != ">": # not a FASTA block, so add this line to the sequence
            entry.seq += line.replace("\n", "").replace("\r", "")

        if line[:1] == ">": # fasta start block
            # start recording
            # add the old Node to the list
            if entry.name != "empty":
                # convert the list to a tuple
                n = entry.seq.lower()
                elementList.append(n) # and stick it on the big list
            entry = Node(line) # make a new node

    # all done;
    return (elementList)

def repeat_mask(seq):
    return(seq.replace("a", "n").replace("c", "n").replace("g", "n").replace("t", "n"))

def loadTSVAsLIST(file):
    oh = open(file, "rU")
    reader = csv.reader(oh, dialect=csv.excel_tab)
    newl = []
    for line in reader:
        newl.append(line)
    return(newl)

def renameDuplicatesFromCSV(path, csvfile, outfile, column_no = 3, bKeepEmptyCols=False):
    """
    append _1 .. _n to duplicates based on the column no
    """
    inf = open(os.path.join(path, csvfile), "rb")
    outf = open(os.path.join(path, outfile), "wb")

    reader = csv.reader(inf)
    writer = csv.writer(outf)
    ulist = []
    nFound = {}

    for line in reader:
        if line[column_no] in ulist:
            # don't write this enty,
            print "Duplicate:", line[column_no]
            if bKeepEmptyCols:
                if column_no == 0:
                    writer.writerow(["%s_%s" % (line[column_no], nFound[line[column_no]])] + line[column_no+1:])
                elif column_no == len(line):
                    writer.writerow(line[:column_no] + ["%s_1" % line[column_no]])
                else:
                    writer.writerow(line[:column_no] + ["%s_1" % line[column_no]] + line[column_no+1:])
                nFound[line[column_no]] += 1
        else:
            # add to ulist and write to file;
            if line[column_no]: # if column is empty don't add it to the list, but write to file
                ulist.append(line[column_no])
                # add a key to the nTimeFound dict;
                nFound[line[column_no]] = 1
                writer.writerow(line) # if I tab this in - don't keep empty rows

            else: # col is empty;
                #print "Duplicate: <empty>"
                if bKeepEmptyCols: writer.writerow(line)
    inf.close()
    outf.close()

"""
    It's hard to believe, but these custom routines below are 10x as
    fast as their numpy equivalents...
    Numpy is a dog.
"""
def mean(intList):
    try:
        return sum(intList) / float(len(intList))
    except TypeError:
        return(intList) # intList is probably a single int -

def std(intList):
    return(math.sqrt(mean([(abs(x - mean(intList) ) ** 2) for x in intList])))

def transpose(list):
    """
    a transpose command, rotates a mtrix or equivalent by 90

    more like a transpose from R than anything else.
    """
    newl = []
    try:
        rows = len(list[0])
    except:
        rows = 1 # probably.
    cols = len(list)

    for row in xrange(rows):
        newl.append([0 for x in xrange(cols)])
    for r in xrange(rows):
        for c in xrange(cols):
            newl[r][c] = list[c][r]
    return(newl)

def isPalindromic(seq):
    """
    is a sequence palindromic?
    returns True or False
    """
    if rc_expanded(seq.lower()) == seq.lower():
        return(True)
    return(False)
