"""
An inherited version of geneList that has extra methods for
dealing with seqfile results.

* this should inherit from peaklist?
"""

import config
from genelist import genelist
from peaklist import peaklist
from draw import draw
from copy import deepcopy

class taglist(genelist):
    def __init__(self, bins, name=""):
        genelist.__init__(self)
        self.draw = draw(self)
        self.name = name
        self.bins = bins

    def saveCSV(self, **kargs):
        print "Error: Saving as a CSV is not available for a 'taglist' (the format is too complex)"
        print "       use the taglist.save() method for a binary save"

    def loadCSV(self, **kargs):
        print "Error: Loading a CSV 'taglist' is not possible, format is too complex."
        print "       Load a taglist using the taglist = load(file) method."

    def _loadCSV(self, **kargs):
        """ Override this internal method """
        pass

    def drawHeatmap(self, **kargs):
        """
        **Purpose**
            draw a heatmap of the taglist.

        **Arguments**

        filename
            filename to save. (full path)

        see also draw._heatmap for other valid arguments.

        **Result**

        * returns True if succesful, None if it fails.
        * Saves a png file to 'filename'
        """
        # arrange the data:
        set_names = []
        col_names = self.bins # comes from seqfile, used to bin the tags.
        row_names = self["name"]
        data = []

        # get colnames
        for name in self["chip_tags"][0]:
            col_names.append(name)

        m = []
        for item in self:
            row = item["chip_tags"]
            for name in row:
                data.append([int(x) for x in row[name]])
                m.append(max(row[name]))
        true_max = max(m)
        #print data

        # add/override any options here.
        kargs["breaks"] = r.seq(0, true_max, float(true_max)/256)
        kargs["data"] = data
        kargs["set_names"] = set_names # why? heatmap.2 does not support this one.
        kargs["col_names"] = col_names
        kargs["row_names"] = row_names
        kargs["colv"] = False

        # pass it up to draw.heatmap()
        self.draw._heatmap(**kargs) # pass the kargs back
        return(True)

    def __add__(self, gene_list):
        """
        (Override)
        This will now also add together the "chip_tag" data.
        """
        # I should test for equivalenvce here.
        newl = self.__copy__()
        newl.linearData = []

        for item in self:
            for other in gene_list:
                if item["name"] == other["name"]:
                    newitem = copy.deepcopy(item)
                    for k in other["chip_tag"]:
                        newitem["chip_tag"][k] = other["chip_tag"][k]
                    newl.linearData.append(newitem)
                else:
                    print "Error: Incompatability in lists."
        newl._optimiseData()
        return(newl)
