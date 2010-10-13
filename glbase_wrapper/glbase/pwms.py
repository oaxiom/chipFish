"""

pwms.py

An amalgamation of a set of pwms.

"""

import sys, os, csv, random
import config, utils

from numpy import zeros, array, mean, std
from draw import draw
from genelist import genelist
from errors import AssertionError, NotImplementedError
from progress import progressbar
from location import location
from history import historyContainer

import pwm

class pwms(genelist):
    def __init__(self, filename=None, name="None", format=None, **kargs):
        """
        **Purpose**
            A collection of pwm items with a few specialised
            methods for the analysis of multiple motifs.

        **Arguments**
            filename (Required)
                some sort of file or path that I can make sense of as a pwm.

            name (Optional)
                The name of the collection of pwms.

            format (Required)
                The format of the pwm file.
                Valid formats currently supported are:
                "uniprobe"
                "transfac"

                to add, but currently unsupported formats:
                "jaspar"
        """
        assert filename, "filename argument not specified"
        assert format, "no format specifier provided"
        assert os.path.exists(filename), "file '%s' not found" % filename

        self.linearData = []
        self.name = name
        self.draw = draw(self)
        self._history = historyContainer(self) # a list of historyItem's

        if format == "uniprobe":
            self.__load_uniprobe(filename)
        elif format == "jaspar":
            raise NotImplementedError
        elif format == "transfac":
            self.__load_transfac(filename)

        self._optimiseData()

    def __load_uniprobe(self, path):
        """
        Load the uniprobe pwm data set.

        To get the uniprobe data set, go to uniprobe>Downloads and
        download 'all data'. You should have a file called "all_data.zip"

        unzip that to give a folder:
        "ALL_DATA/"

        this function will then step through the data and load all of it in.
        """
        for subdir, dirs, files in os.walk(path):
            for file in files:
                # get rid of the extension and use this to append the motif.
                if (file[0] != ".") and ("readme" not in file): # Mac OSX specific malarky.
                    new_pwm = {"name": file.split(".")[0]}

                    oh = open(os.path.join(subdir, file), "rU")

                    name = oh.readline()
                    data = []
                    for line in oh:
                        l = line.split("\t")[1:]
                        l = [float(i) for i in l]
                        if len(l) > 0:
                            data.append(l)

                    # convert to numpy.
                    pfm = array(data)
                    pfm = pfm.T

                    self.linearData.append({"name": new_pwm["name"], "pwm": pwm.pwm(name=new_pwm["name"], pwm_matrix=pfm, isPFM=True)})

    def __load_transfac(self, filename):
        """
        Load in transfac db. In this case transfac comes as a single file,
        so you need to send me that file.
        """
        oh = open(filename, "rU")

        newpwm = {"desc":""} # temp store.
        data = []
        degen = []

        for line in oh:
            two_chars = line[0:2]

            # process easily read labels:
            if two_chars == "ID":
                newpwm["ID"] = line[2:].replace("\n", "")
            elif two_chars == "NA":
                newpwm["name"] = line[2:].replace("\n", "")
            elif two_chars == "DE":
                newpwm["desc"] = line[2:].replace("\n", "")
            elif two_chars == "//":
                # finish the pwm, clean it up and generate a new one.
                pfm = array(data)
                self.linearData.append({"id": newpwm["ID"], "pwm": pwm.pwm(name=newpwm["name"], pwm_matrix=pfm, isPFM=False),
                    "degen": "".join(degen), "name": newpwm["name"], "desc": newpwm["desc"]}) # don't convert to pwm
                newpwm = {"desc": ""}
                data = []
                degen = []

            # try to grab the motif:
            try:
                p = int(line[:2])
                while p > len(data):
                    data.append([])
                while p > len(degen):
                    degen.append("")

                t = line.replace("\n", "").split(" ")[1:]
                newt = [float(i) for i in t[:-1] if i != ""]

                data[p-1] = newt
                degen[p-1] = t[-1]
            except ValueError:
                pass

        return(None)

    def _optimiseData(self):
        """
        (Internal)
        not used here.
        """
        pass

    def __repr__(self):
        return("glbase.pwms")

    def scan_sequence(self, loc, sequence, features, filename=None,
        merge_strands=True, summary_file=None,
        z_score=4.5, **kargs):
        """
        **Purpose**

        **Arguments**
            sequence
                The DNA sequence to scan.

            features
                send a dictionary, derived from a genelist loaded from
                say refseq/ucsc and this will know what to do with it.

            filename (Required)
                The filename to save the image files to.

            summary_filename (Required)
                filename to save a summary of the interesting hits to,
                saves a tsv file.

            merge_strands (Optional, default=True)
                merge the score strands. At the moment only this method
                is supported and separate scores per strand are not supported.

            z_score (Optional, default = 4.0)
                The z-score to use to call an interesting peak.
                By default it is set to 4 standard deviations
                away from the mean.

        **Returns**
            saves an image file to filename.
        """
        result = []

        # cooerce the location
        loc = location(loc=loc)

        p = progressbar(len(self))
        for pi, pwm in enumerate(self.linearData):
            data = pwm["pwm"].scan_sequence(sequence, merge_strands)

            # get the significant hits
            m = mean(data)
            s = std(data)
            # collect interesting results
            z_reqd = m+(s*z_score)
            for i, v in enumerate(data):
                if v > z_reqd:
                    result.append({"loc": location(chr=loc["chr"], left=loc["left"]+i, right=loc["left"]+i+len(pwm["pwm"])),
                        "score": v, "strand": "?", "seq": sequence[i:i+len(pwm["pwm"])], "Z-score": abs(m-v)/s,
                        "name": pwm["name"], "__local_location": i})

            real_filename = "%s_%s.png" % (self.name, pwm["name"])
            self.draw._plot_and_histogram(filename=real_filename, data=data,
                figsize=(20,5), title=pwm["name"], ylims=(0,1), xlims=(0, len(data)),
                loc=loc, genomic_features=features)

            del data # clean up the array, some really big searches kill
                     # the 'puter.

            p.update(pi)
        config.log.info("Saved %s image files" % len(self))

        oh = open(summary_file, "w")
        oh.write("motif\tloc\tscore\tstrand\tseq\tZ-score\n")
        for item in result:
            oh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (item["name"], str(item["loc"]), item["score"],
                item["strand"], item["seq"], item["Z-score"]))
        oh.close()
        config.log.info("Saved a summary file to '%s'" % summary_file)

        # save a summary label image.
        sum = []
        for item in result:
            sum.append({"pos": (item["__local_location"], item["Z-score"]), "label": item["name"]})
        real_filename = self.draw._labeled_figure(data=sum, axissize=(len(sequence), 0), ylim=(4,5.5),
            figsize=(20,5), genomic_features=features, filename=summary_file,
            loc=loc)
        config.log.info("Saved a summary image file to '%s'" % real_filename)

        return(real_filename)

    def saveCSV(self, filename=None, **kargs):
        """
        **Purpose**
            save the pwm list as a csv.
            This doens't work completely - the pwm matrix is not correctly saved.

        **Arguments**
            see genelist.saveCSV()
        """
        if "key_order" in kargs and kargs["key_order"]:
            key_order = kargs["key_order"]
            if "pwm" in key_order:
                key_order.remove("pwm")
        else:
            key_order = [k for k in self.linearData[0] if k != "pwm"]
        key_order.append("pwm")

        genelist.saveCSV(self, filename=filename, key_order=key_order, **kargs) # inherit
