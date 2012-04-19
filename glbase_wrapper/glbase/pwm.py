"""

pwm.py

Tools and utilities to use position-weight matrices and other matrix-like
representations

TODO:
-----
. Merge with logo.py

"""

from __future__ import division

import config

import sys, os, math
from numpy import zeros, array

from base_genelist import _base_genelist
from genelist import genelist
from draw import draw
from utils import rc, convertFASTAtoDict
from progress import progressbar
from location import location

# see if weblogo is available.
try:
    import weblogo
    WEBLOGO_AVAILABLE = True
except:
    WEBLOGO_AVAILABLE = False # fail silently

class pwm:
    """
    **Purpose**
        Store a definition of a pwm_matrix

    **Arguments**
        name (Required)
            name of the matrix. It is required because many downstream
            scanning efforts require a unique name.

        pwm_matrix (Required, or fasta_file)
            a matrix describing the pwm. For example:
            [ [a, c, g, t],
              [a, c, g, t],
              ...
              [a, c, g, t] ]
            or you can use a numpy 2D array as well.

            You can also send the matrix as a list of dictionaries.
            [ {"a": 1, "c": 2, "g": g, "t": t} ... ]

        fasta_file (Required, or pwm_matrix)
            load a fasta file and turn int into a pwm_matrix
            
        txt_file (Required)
            load a txt file which is a list of sequences and convert into a 
            pwm/pfm.

        isPFM (Optional, default=True)
            the matrix is actually a frequency matrix that needs to be
            converted into a weight matrix. Set to False if your pwm is already a log-odds 
            position weight matrix.
    """
    def __init__(self, name, pwm_matrix=None, fasta_file=None, txt_file=None, isPFM=True):

        self.name = name
        
        # get around numpy.array() __non_zero__ silliness:
        try:
            __load_matrix = False
            if pwm_matrix:
                __load_matrix = True
        except ValueError: 
            if pwm_matrix.any():
                __load_matrix = True
        
        if __load_matrix:
            if isinstance(pwm_matrix, list):
                # convert to a numpy array
                if isinstance(pwm_matrix[0], dict):
                    self.__matrix = self.__convert_dict_matrix(self, pwm_matrix)
                else:
                    l = len(pwm_matrix[0])
                    pfm = array( pwm_matrix , dtype=float)
                    self.__matrix = pfm
            else:
                self.__matrix = pwm_matrix # Probably a numpy array?
        elif fasta_file:
            g = convertFASTAtoDict(fasta_file)
            pwm_matrix = None
            for e in g:
                if not pwm_matrix: # matrix sizes not sampled yet
                    s = len(e["seq"])
                    pwm_matrix = [{"a": 0, "c": 0, "g": 0, "t": 0} for n in range(s)]
                
                for i, bp in enumerate(e["seq"]):
                    pwm_matrix[i][bp] += 1
            self.__matrix = self.__convert_dict_matrix(pwm_matrix)
        elif txt_file:
            oh = open(txt_file, "rU")
            for line in oh:
                if not ">" in line:
                    e = line.lower().strip()
                    if not pwm_matrix: # matrix sizes not sampled yet
                        s = len(e)
                        pwm_matrix = [{"a": 0, "c": 0, "g": 0, "t": 0} for n in range(s)]
                    
                    for i, bp in enumerate(e):
                        pwm_matrix[i][bp] += 1
            self.__matrix = self.__convert_dict_matrix(pwm_matrix)
        else:
            raise AssertionError, "no valid input found for pwm"
            
        if isPFM:
            config.log.debug("Converting the frequency matrix '%s' to a log-odds weight matrix" % self.name)
            self.convertPFMtoPWM() 
        self.__do_minmax()

    def __convert_dict_matrix(self, dict_pwm_matrix):
        """ 
        load in the matix when it comes in the form: [ {"a": 1, "c": 2, "g": g, "t": t} ... ]
        
        Convert to:                
            [[a, c, g, t],
            [a, c, g, t],
            ...
            [a, c, g, t] ]
        """
        new_matrix = []
        for i in dict_pwm_matrix:
            new_matrix.append([i["a"], i["c"], i["g"], i["t"]])
        return(array(new_matrix , dtype=float))

    def convertPFMtoPWM(self):
        """
        **Purpose**
            Convert a frequency matrix into a weight matrix.

            Converts a pfm to a pwm
            matrices must be in pwm format for accurate searching.
            via this algorithm:

            w = log2 ( ( f + sqrt(N) * p ) / ( N + sqrt(N) ) / p )

            where:
            w - is a weight for the current nucleotide we are calculating
            f - is a number of occurences of the current nucleotide in the current column (e.g., "1" for A in column 1, "8" for C etc)
            N - total number of observations, the sum of all nucleotides occurences in a column (13 in this example)
            p - [prior] [background] frequency of the current nucleotide; this one usually defaults to 0.25 (i.e. one nucleotide out of four)

            Shamelessly borrowed from:
            http://bogdan.org.ua/2006/09/11/position-frequency-matrix-to-position-weight-matrix-pfm2pwm.html

        **Arguments**
            None

        **Results**
            This object now contains a weight matrix, not a frequency matrix.
        """
        for row in xrange(len(self.__matrix)):
            bign = sum(self.__matrix[row])
            for bp in xrange(len(self.__matrix[row])):
                self.__matrix[row][bp] = math.log( ( self.__matrix[row][bp] + math.sqrt(bign) * 0.25) / (bign + math.sqrt(bign)) / 0.25, 2) # log 2, not log n

        self.__do_minmax() # redo min/max scores

    def __len__(self):
        """
        Method to give a valid len(self) assignment.
        """
        return(len(self.__matrix)) # just pass back the length

    def __do_minmax(self):
        """
        Redo the min max score (for example, after doing pfm -> pwm)
        """
        n = 0
        self.__minscore = 0
        self.__maxscore = 0
        for n in self.__matrix:
            self.__minscore += min(n)
            self.__maxscore += max(n)

    def __str__(self):
        return("<name: %s length: %s, minmax: (%.1f, %.1f)>" % (self.name, self.__len__(), self.__minscore, self.__maxscore))

    def get_matrix(self):
        """
        **Purpose**
            Return the actual pwm matrix.
            THis is a numpy array and is the native format for the matrix
        """
        return(self.__matrix)

    def get_as_list(self):
        """
        **Purpose**
            Return the pwm as a list, in the form [ [a], [c], [g], [t] ]
            (This can be loaded straight into MOODS
        """
        c = self.__matrix.T # get a transposed version.
        new = []
        new.append([v for v in c[0]]) # a # yes, yes, could be done on a single line...
        new.append([v for v in c[1]]) # c # for clarity.
        new.append([v for v in c[2]]) # g
        new.append([v for v in c[3]]) # t
        return(new)

    def score(self, sequence):
        """
        **Purpose**
            return the pwm threshold score for a particular sequence.

        **Arguments**
            sequence
                a dna sequence string
                Note that this method only scans the sequence for the
                length of the pwm. It will throw an exception if the
                length of the sequence is different than the length of the motif.
                (This behaviour is intentional)

                If the seq contains any N nucleotides the score will be
                returned as 0.0

        **Returns**
            A dictionary {"+": <upper strand score>, "-": <lower strand score>}
        """
        seq = sequence.lower()

        if seq.count("n") > 0:
            return({"+": 0.0, "-": 0.0}) # reject if contains an N

        con_setpos = {"a" : 0, "c": 1, "g": 2, "t": 3} # convert dict location to matrix location
        result = {} # new list

        seq_data = {"+": seq, "-": rc(seq)}
        for key in seq_data: # new super small version:
            score_list = [self.__matrix[i][con_setpos[letter]] for i, letter in enumerate(seq_data[key])]
            unnormalised_score = sum(score_list)
            result[key] = (unnormalised_score - self.__minscore) / (self.__maxscore - self.__minscore)

        return(result)

    def scan_sequence(self, sequence, merge_strands=True):
        """
        **Purpose**
            Scan across a sequence and return the score as an array.

        **Arguments**
            sequence (Required)
                The sequence to search. A string

            merge_strands (Optional, default=True)
                Each strand provides a seperate binding site.
                This can be reported seperately or the
                strands can be merged, reporting only the highest
                threshold score for that site.

                If True, returns a single list.
                If False returns a dictionary {"+": [0...n], "-": [0...n]}

        **Returns**
            By default a list, but see Arguments, merge_strands for
            more details.
        """
        if merge_strands:
            result = zeros(len(sequence)-len(self))
        else:
            result = {"+": zeros(len(sequence)), "-": zeros(len(sequence))}

        for p in xrange(len(sequence)-len(self)):
            scores = self.score(sequence[p:p+len(self)])

            if merge_strands:
                result[p] = max([scores["+"], scores["-"]])
            else:
                result["+"][p] = scores["+"]
                result["-"][p] = scores["-"]

        return(result)

    def scan_seq_with_features(self, sequence, features=None, filename=None, **kargs):
        """
        **Purpose**
            scan a sequence and produce a graph with annotated features.

        **Arguments**
            sequence
                A string containing a sequence of a dict/genelist item
                with a "loc" key

            features
                a dictionary containing a location and a label.

            filename
                filename to save the image to.

        **Returns**
            The actual filename used to save the figure.
        """
        raise NotImplementedError

    def scan_genelist(self, fasta_list, threshold=0.75, keep_empty=True):
        """
        **Purpose**
            Scan a fasta list (i.e. a genelist with a 'seq' key and
            corresponding sequence data)

        **Arguments**
            fasta_list
                A fasta list.
            
            threshold
                The threshold to pass the motif on
             
            keep_empty (Optional, default=False)
                Keep fasta entries that do not have a motif

        **Returns**
            A new genelist-like object, with these new key:
                "pos", "strand", "seq", "score", "motifs" and the "name" from the original FASTA list
            If the original list also has a "loc" key containing a valid genomic location then
            it will contain the genomic location of the motif stored in the new key "motif_loc"
        """
        newl = genelist()
        newl.name = fasta_list.name
                
        p = progressbar(len(fasta_list))
        for n, item in enumerate(fasta_list):
            res = self.scan_sequence(item["seq"], False)
           
            if res:
                seq_data = {"+": item["seq"], "-": rc(item["seq"])}
                for strand in ("+", "-"):
                    for i, score in enumerate(res[strand]):
                        if score > threshold:
                            newi = {}
                            newi["name"] = item["name"]
                            newi["local_pos"] = i
                            if "loc" in item:
                                if strand == "+":
                                    newi["motif_loc"] = location(chr=item["loc"]["chr"], left=item["loc"]["left"]+i, right=item["loc"]["left"]+i+len(self))
                                else:
                                    newi["motif_loc"] = location(chr=item["loc"]["chr"], left=item["loc"]["right"]-i-len(self), right=item["loc"]["right"]-i)
                                
                            newi["strand"] = strand
                            if strand == "+":
                                newi["seq"] = seq_data[strand][i:i+len(self)]
                            else:
                                newi["seq"] = rc(seq_data[strand][i:i+len(self)])
                            newi["score"] = score
                            newi["motifs"] = "Yes"
                            newl.linearData.append(newi)
            else:
                newl["motifs"] = None
            p.update(n)
        
        newl._optimiseData()
        return(newl)

    def save(self, filename=None, mode="homer"):
        """
        **Purpsose**
            save the pwm into another file.
            
        **Arguments**
            filename (Required)
                The filename to save the motif to.
        
            mode (Optional, default="homer")
                the type or mode of saving the file. 
                At the moment these formats are supported:
                
                homer - save as a pwm
        
        **Returns**
            None
        """
        assert filename, "no filename specified"
        
        if mode == "homer":
            oh = open(filename, "w")
            
            oh.write(">%s\t%s\t%s\t%s\t%s\t%s\n" % (self.name, self.name, 0, 0, 0, "T:0(0),B:0(0),P(0)"))
            for i in self.__matrix:
                nl = i/sum(i)
                oh.write("%s\n" % "\t".join([str(b) for b in nl]))          
        elif mode == "counts":
            oh = open(filename, "w")
            
            oh.write(">%s\t%s\t%s\t%s\t%s\t%s\n" % (self.name, self.name, 0, 0, 0, "T:0(0),B:0(0),P(0)"))
            for i in self.__matrix:
                oh.write("%s\n" % "\t".join([str(b) for b in nl])) 
        
        return(None)

if __name__ == "__main__":
    # These need to be moved into tests:

    m = pwm("soxoct",
            [[4,    43, 0,  8],         # c -sox2 start
            [37,    0,  0,  25],        # a/t
            [0,     0,  0,  79],        # t
            [4,     11, 0,  64],        # t
            [14,    2,  54, 8],         # g
            [9,     3,  17, 62],        # t -sox2 end
            [13,    58, 4,  24],        # c/g #
            [99,    0,  0,  0],         # a - oct4 start
            [0,     0,  0,  99],        # t
            [4,     0,  84, 12],        # g
            [42,    42, 4,  4],         # a/c
            [40,    39, 0,  4],         # a/c
            [56,    3,  12, 2],         # a
            [62,    0,  0,  1],         # a
            [16,    0,  2,  29]])
    print m
    print m.score("cattgtcatgcaaat")
    print m.score("cattgacacgcttat")

    print m.scan_sequence("aaaaaaaaaatttgcatgacaagaagccttcttttttttttcttcattgtcatgcaaatcccccggggg") # two sox-oct motifs.

    print
    print m.get_matrix()

    m = pwm("meh", txt_file="tests/txt_motif_file.txt", isPFM=False)

    m.save("meh.matrix")