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
            res.append(re.compile(newSeq))
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

        Scans both strands, but returns only a single score per sequence
        in the list (ie score = forward + reverse).

        **Arguments**

        genelist
            genelist with sequence data attached, must be in a "seq", "sequence"
            "f". If you load a fasta, the sequence is loaded into "seq".
            You can send a list of genelists as well, and these will all be plotted
            on the graph.

        random_fasta (Optional, default=None)
            a random or list of random_fasta lists to act as a comparison.

        filename
            filename of the resulting image file.

        window (Optional, default = 10% of the first list in 'genelist')
            the moving window average size

        resample_random (Optional, default = True)
            resample the random fasta lists so that they equal the size
            of the original list. Set to False if you don't want to do this.

        **Result**

        returns None and a figure saved to filename
        """
        valid_args = ["genelist", "random_fasta", "filename", "window", "resample_random"]
        for key in kargs:
            assert key in valid_args, "scanMotifFrequency() - Argument '%s' not recognised" % key

        # draw a moving average plot
        assert filename, "you must specify a filename to save the graph to"
        assert genelist, "you must provide a genelist"

        # turn genelist into a vanilla list if it is not already.
        if not isinstance(genelist, list):
            genelist = [genelist]

        window = int(len(genelist[0]) * 0.10) # 10% of the list
        if kargs.has_key("window"):
            window = int(kargs["window"])

        motif = self.getRegEx()
        __bWarningDoneAlready = False

        if kargs.has_key("random_fasta"):
            random_fastas = kargs["random_fasta"]
            if kargs.has_key("resample_random") and (not kargs["resample_random"]):
                pass # resample_random = True by default
            else: # default resample:
                new_fasta_list = []
                for item in random_fastas:
                    if len(genelist[0]) > len(item):
                        # just print a warning
                        # Although, make sure only to print it once so I don't look stupid.
                        if not __bWarningDoneAlready:
                            print "Warning: The genelist is larger than the random lists, resampling will result in duplidate entries"
                            __bWarningDoneAlready = True
                        list_size = len(item)
                    else:
                        list_size = len(genelist[0])

                    # and then resample anyway
                    newl = new_gl() # get a new and empty list # newgl is an empty genelist

                    while len(newl) != len(genelist[0]): # don't use random.sample() as this can bodge larger sample lists (naughty!!)
                        r = random.randint(0, list_size-1)
                        newl.linearData.append(item.linearData[r])
                    newl._optimiseData()
                    new_fasta_list.append(newl)

                random_fastas = new_fasta_list
            fasta_lists = genelist + random_fastas
        else:
            fasta_lists = genelist

        # get the seq key:
        if "seq" in genelist[0].getKeys():
            seq_key = "seq"
        elif "f" in genelist[0].getKeys():
            seq_key = "f"
        elif "sequence" in genelist[0].getKeys():
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

        actual_filename = self.draw._qplotxy(motif_result, filename=filename, labels=labels)
        print "Info: Saved figure %s" % actual_filename
