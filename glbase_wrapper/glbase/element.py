"""

element.py

lists and tools to deal with motifs.

Part of glbase.

"""

import config

import re, random, sys, os, string

import utils
from genelist import genelist as new_gl
from errors import AssertionError
from draw import draw

# R=[AG], Y=[CT], K=[GT], M=[AC], S=[GC], W=[AT],
# [ACT] = h, [ACG] = v, [AGT] = d, [CGT] = b
regex_dict = {
    "a" : "a",
    "c" : "c",
    "g" : "g",
    "t" : "t",
    "r" : "[ag]", # twos
    "y" : "[ct]",
    "k" : "[gt]",
    "m" : "[ac]",
    "s" : "[gc]",
    "w" : "[at]",
    "h" : "[act]", # threes
    "v" : "[acg]",
    "d" : "[agt]",
    "b" : "[cgt]",
    "n" : "[acgt]" # four
}

class motif:
    def __init__(self, name=None, sequence="", **kargs):
        """
        constructor for an element
        an element here is a piece of DNA proposed to be a TF binding site.
        for example:
        you can make an empty element with:
        element("empty_name")
        eg:
        Sox2: ywttswn
        Oct4: atgywnww
        """
        assert sequence, "no specified sequence for the motif"
        assert name, "you must give your motif a name"

        self.name = name
        self.draw = draw(self)
        self.seq = sequence.lower()
        self.palindromic = utils.isPalindromic(sequence)

        if kargs.has_key("scramble") and kargs["scramble"]: self._scrambleMotif()

    def isPalindromic(self):
        return(self.bPalindromic)

    def getRegEx(self):
        """
        return the regex.
        """
        self.re = re.compile("".join([regex_dict[bp] for bp in self.seq]))
        return(self.re)

    def _scrambleMotif(self, number=10):
        """
        (Internal)
        generate scrambled motifs.

        Be careful with this, it does not work the same as in fexcom.
        It does not scramble the motif in-place, instead it returns a set of
        'regexd' motifs.
        """
        res = []
        for n in xrange(number):
            newSeq = oldSeq = self.seq
            while newSeq == oldSeq: # stops the new motif == old motif.
                q = []
                for item in self.seq:
                    q.append(item)

                new = []
                while q:
                    i = random.randint(0, len(q)-1)
                    new.append(q[i])
                    q.remove(q[i])

                newSeq = "".join(new)
                print new
            res.append(re.compile(newSeq))
        print res
        return(res)

    def __str__(self):
        return("motif name:%s sequence:%s" % (self.name, self.seq))

    def __len__(self):
        return(len(self.seq))

    def scanMotifFrequency(self, genelist=None, filename=None, **kargs):
        """
        **Purpose**

        scan a set of sequences for the maximum score for the current motif.
        score each fasta sequence for a match for a particular motif
        then output the result as a movingAverage graph.

        **Arguments**

        genelist
            genelist with sequence data attached, must be in a "seq", "sequence"
            "f", or similar key.

        random_fasta (Optional, default=None)
            a random or list of random_fasta lists to act as a comparison.

        filename
            filename of the resulting image file.

        window (Optional, default = 5%)
            the moving window average size

        resample_random (Optional, default = True)
            resample the random fasta lists so that they equal the size
            of the original list.

        **Result**

        returns None and a figure saved to filename
        """
        # draw a moving average plot
        assert filename, "you must specify a filename to save the graph to"
        assert genelist, "you must provide a genelist"

        window = int(len(genelist) * 0.10) # 5% of the list
        if kargs.has_key("window"):
            window = int(kargs["window"])

        motif = self.getRegEx()
        bWarningDoneAlready = False

        if kargs.has_key("random_fasta"):
            fasta_lists = kargs["random_fasta"]
            if kargs.has_key("resample_random") and (not kargs["resample_random"]):
                pass # resample_random = True by default
            else: # default resample:
                new_fasta_list = []
                for item in fasta_lists:
                    if len(genelist) > len(item):
                        # just print a warning
                        # Although, make sure only to print it once so I don't look stupid.
                        if not bWarningDoneAlready:
                            print "Warning: The genelist is larger than the random lists, resampling will result in duplidate entries"
                            bWarningDoneAlready = True
                        list_size = len(item)
                    else:
                        list_size = len(genelist)

                    # and then resample anyway
                    newl = new_gl() # get a new and empty list # newgl is an empty genelist

                    while len(newl) != len(genelist):
                        r = random.randint(0, list_size-1)
                        newl.linearData.append(item.linearData[r])
                    newl._optimiseData()
                    new_fasta_list.append(newl)

                fasta_lists = new_fasta_list
            fasta_lists = [genelist] + fasta_lists
        else:
            fasta_lists = [genelist]

        # get the seq key:
        if "seq" in genelist.getKeys():
            seq_key = "seq"
        elif "f" in genelist.getKeys():
            seq_key = "f"
        elif "sequence" in genelist.getKeys():
            seq_key = "sequence"
        else:
            raise AssertionError, "the genelist does not appear to have a sequence attached" # fake an assert

        motif_result = []
        labels = []

        for fasta in fasta_lists:
            res = []
            for item in fasta.linearData:
                matches = []
                for seq in [item[seq_key], utils.rc(item[seq_key])]:
                    matches += motif.findall(seq.lower())
                res.append(len(matches))
            labels.append(fasta.name)
            motif_result.append(utils.movingAverage(res, window, False))

        actual_filename = self.draw._qplot(motif_result, filename=filename, labels=labels)
        print "Info: Saved figure %s" % actual_filename
