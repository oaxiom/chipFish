"""

glglobs are almagamated lists, useful for drawing comparisons between lists.

glglobs also inherit from genelist, but many of the methods are broken.

renamed to glglob as it clashes with a matplotlib and python module name

"""

from __future__ import division

import sys, os, csv, string, math

from numpy import array, zeros, object_, arange
from array import array as qarray
from copy import deepcopy

import config, utils
from flags import *
from genelist import genelist
from draw import draw
from errors import AssertionError, NotImplementedError, GlglobDuplicateNameError

import matplotlib.pyplot as plot
import matplotlib.cm as cm
from scipy.stats import spearmanr, pearsonr

class glglob(genelist): # cannot be a genelist, as it has no keys...
    def __init__(self, *args, **kargs):
        genelist.__init__(self)

        # args should be a list of lists.
        # we then store them in the linearData set
        #for a in args: # no need for this restriction?
        #    for b in args:
        #        assert a == b, "Lists must be equivalent"
        self.linearData = args
        self._optimiseData()

    def __repr__(self):
        return("glbase.glglob")

    def __str__(self):
        # work out all of the types
        types = []
        for item in self.linearData:
            if item.__repr__() not in types:
                types.append(item.__repr__())
        return("glglob contains: %s items of type(s): %s" % (len(self.linearData), ", ".join(types)))

    def _optimiseData(self): # no keys, so would die.
        """
        (Override)

        actually does have keys: the list names
        """
        self.__list_name_lookback = {} # the index location.
        for index, item in enumerate(self.linearData):
            if item.name in self.__list_name_lookback:
                raise GlglobDuplicateNameError, ("self._optimiseData", item.name)
            else:
                self.__list_name_lookback[item.name] = index

    def loadCSV(self):
        config.log.error("glglobs cannot be sensibly represented as a CSV file, use glload() to load binary copies")
        return(False)

    def saveCSV(self):
        config.log.error("glglobs cannot be sensibly represented as a CSV file, use .save() to save binary copies")
        return(False)

    def __getitem__(self, value):
        """
        (Override)

        glglobs should be accesible by name.
        """
        if value in self.__list_name_lookback:
            return(self.linearData[self.__list_name_lookback[value]]) # get the actual item
        return(None)

    def __setitem__(self, value):
        """
        (Override)
        """
        config.log.error("glglobs cannot be written to")

    def compare(self, key=None, filename="", method=None, **kargs):
        """
        **Purpose**
            perform a square comparison between the genelists according
            to some sort of criteria.

        **Arguments**
            key (string, required)

            filename (string, required)
                filename to save heatmap image as.

            method (string, "overlap|collide|map", required)
                method to use to compare the genelists

            distance_score (string, "euclidean", optional defaults to "euclidean")
                Scoring method for distance caluclations.
                This is not implemnted at the moment, and only uses a
                euclidean distance.

            output_pair_wise_correlation_plots (True|False)
                This is a debug option to output graphs showing the
                actual correlation matrices of the pair-wise correlations.
                Probably best only to do this if you
                know what you are doing!

        **Result**
            returns the distance matrix if succesful or False|None if not.
            Saves an image to 'filename' containing a grid heatmap
        """
        valid_methods = ["overlap", "collide", "map"]
        valid_dist_score = ["euclidean"]
        distance_score = "euclidean" # override for the moment
        assert method in valid_methods, "must use a valid method for comparison (%s)" % ", ".join(valid_methods)
        assert filename, "Filename must be valid"
        assert key in self.linearData[0].linearData[0], "key '%s' not found" % (key,) # just check one of the lists
        assert distance_score in valid_dist_score, "%s is not a valid distance metric" % (distance_score,)

        config.log.info("This may take a while, all lists are intersected by '%s' with '%s' key" % (method, key))

        matrix = zeros( (len(self), len(self)) ) # 2D matrix.

        for ia, la in enumerate(self.linearData):
            for ib, lb in enumerate(self.linearData):
                if ia == ib:
                    matrix[ia, ib] = len(la) # should be filled in with the maximum possible overlap.
                elif ia < ib: # make search triangular
                    pass
                else:
                    if method == "collide":
                        matrix[ia, ib] = len(la.collide(genelist=lb, loc_key=key, delta=200))
                    elif method == "overlap":
                        matrix[ia, ib] = len(la.overlap(genelist=lb, loc_key=key, delta=200))
                    elif method == "map":
                        matrix[ia, ib] = la.map(genelist=lb, key=key)

        # fill in the gaps in the triangle
        for ia, la in enumerate(self.linearData):
            for ib, lb in enumerate(self.linearData):
                if ia < ib:
                    matrix[ia,ib] = matrix[ib,ia]

        # data must be normalised to the maximum possible overlap.
        for ia, la in enumerate(self.linearData):
            for ib, lb in enumerate(self.linearData):
                matrix[ia,ib] = (matrix[ia,ib] / min([len(la), len(lb)]))

        #print matrix
        #print spear_result_table

        spear_result_table = zeros( (len(self), len(self)) ) # square matrix to store the data.
        # convert the data to a spearmanr score.
        for ia, this_col in enumerate(matrix):
            for ib, other_col in enumerate(matrix):
                #if ia != ib:
                #spear_result_table[ia,ib] = spearmanr(this_col, other_col)[0] # [0] = r score, [1] = p-value
                spear_result_table[ia,ib] = pearsonr(this_col, other_col)[0] # [0] = r score, [1] = p-value

        #print spear_result_table
        result_table = spear_result_table

        """
        # for euclidean distanes caluclate the distance for each column.
        # Although the data is actually triangluar I fill the whole lot in.
        # This may need optimising later.
        # see the test_glglob for a test example.
        # this is the method for euclidean distance.
        # get comparisons by column.
        # There is probably a much easier way of doing this...
        # slicing the matrices...
        maximum_value = 0 # I want to store this for bracketing the heatmap
        minimum_value = 2147483646 # start with a huge value. = LongInt - 1
        for ci1, this_col in enumerate(matrix):
            for ci2, other_col in enumerate(matrix):
                pair_distance = 0.0
                for each_pair in xrange(len(this_col)):
                    pair_distance += math.pow(float(this_col[each_pair] - other_col[each_pair]), 2)
                result_table[ci1, ci2] = math.sqrt(pair_distance)
                if ci1 != ci2: # ignore the diagonal.
                    if result_table[ci1, ci2] > maximum_value: # I want to store this to bracketing the heatmap
                        maximum_value = result_table[ci1, ci2]
        """
        #result_table = spear_result_table # zeros( (len(self), len(self)) ) # square matrix to store the data.

        # need to add the labels and serialise into a doct of lists.
        dict_of_lists = {}
        row_names = []
        for index, item in enumerate(self.linearData):
            dict_of_lists[item.name] = result_table[index]
            row_names.append(item.name) # preserve order of row names.

        if "output_pair_wise_correlation_plots" in kargs and kargs["output_pair_wise_correlation_plots"]:
            # output the plot matrices.
            for ia, la in enumerate(self.linearData):
                for ib, lb in enumerate(self.linearData):
                    plot_data = []
                    if ia != ib:
                        x = matrix[ia,]
                        y = matrix[ib,]
                        self.draw._scatter(x, y, xlabel=row_names[ia], ylabel=row_names[ib],
                            filename="dpwc_plot_%s_%s.png" % (row_names[ia], row_names[ib]))

        # draw the heatmap and save:
        realfilename = self.draw._heatmap(data=dict_of_lists, filename=filename,
            colbar_label="correlation", bracket=[-0.2, 1],
            square=True, cmap=cm.hot, cluster_mode="euclidean",
            row_names=row_names, col_names=row_names)

        config.log.info("Saved Figure to '%s'" % realfilename)
        return(result_table)

    def venn_diagrams(self, key=None, filename=None, **kargs):
        """
        **Purpose**
            draw 2,3,4 venn Diagrams
            currently only equally size venndiagrams are supported.
            (proportional venn diagrams are experimental only, enable them
            using experimental_proportional_venn = True as an argument).

            your glglob should be loaded with several genelist-like objects
            that have a valid map() function (i.e. pretty much all of them)

            Note that you can do simple 2 overlap venn_diagrams using any
            pair of genelists with this sort of code:

            genelist.map(genelist=other_genelist, <...>, image_filename="venndiagram.png")

        **Arguments**
            key
                key to use to map() between the two lists.

            filename
                save the venn diagram to this filename.

            title (Optional)
                title for the figures
                defaults to <list> vs <list> vs ...

        **Returns**
            A venn diagram saved in filename.
        """
        valid_args = ["filename", "key", "title", "experimental_proportional_venn"]
        for k in kargs:
            if not k in valid_args:
                raise ArgumentError, (self.map, k)

        assert len(self.linearData) <= 4, "currently glglob venn diagrams only support at most 4 overlaps"
        assert len(self.linearData) >= 2, "you must send at least two lists"
        assert key, "Must specify a 'key' to map the two lists"
        assert filename, "no filename specified for venn_diagrams to save to"

        proportional = False
        if "experimental_proportional_venn" in kargs and kargs["experimental_proportional_venn"]:
            proportional=True

        if len(self.linearData) == 2:
            # simple case, genelist.map() can deal with this one.
            self.linearData[0].map(genelist=self.linearData[1], image_filename=kargs["filename"], venn_proportional=proportional)
            config.log.info("Save a venn diagram to: %s" % kargs["filename"])
            return(None)
        elif len(self.linearData) == 3:
            raise NotImplementedError, "Hillariously, 2 and 4 venn diagrams are implemented. But not 3..."
        elif len(self.linearData) == 4:
            config.log.warning("Due to the nature of 4-way venn diagrams, opposing corners are not overlapped as a pair")
            config.log.warning("You should be aware of this limitation")
            # work out the overlap matrix:

            A = self.linearData[0] # A # this is more for clarity than neccesity
            B = self.linearData[1] # B
            C = self.linearData[2] # C
            D = self.linearData[3] # D

            # calculate first round of overlaps (twos):
            AC = A.map(genelist=C, key=key) # AC
            AB = A.map(genelist=B, key=key) # AB
            CD = C.map(genelist=D, key=key) # CD
            BD = B.map(genelist=D, key=key) # BD
            #AD = A.map(genelist=D, key=key) # AD # these two are not presented.
            #BC = B.map(genelist=C, key=key) # BC

            # tidy up any genelists with no entries.
            # second round (threes)
            # fake some lists:
            ABC = []
            ACD = []
            ABD = []
            BCD = []
            if AC and AB:
                ABC = AB.map(genelist=AC, key=key) # ABC
            if AC and CD:
                ACD = AC.map(genelist=CD, key=key) # ACD
            if AB and BD:
                ABD = AB.map(genelist=BD, key=key) # ABD
            if BD and CD:
                BCD = BD.map(genelist=CD, key=key) # BCD

            # FOUR:
            ABCD = []
            if ABC and BCD:
                ABCD = ABC.map(genelist=BCD, key=key) # or any other pair will do.

            lists = {"A": A, "B": B, "C": C, "D": D}
                #"AC": AB, "AB": AB, "CD": CD, "BD": BD,
                #"ABC": ABC, "ACD": ACD, "ABD": ABD, "BCD": BCD,
                #"ABCD": ABCD}

            scores = {"A": len(A), "B": len(B), "C": len(C), "D": len(D),
                "AC": len(AB), "AB": len(AB), "CD": len(CD), "BD": len(BD),
                "ABC": len(ABC), "ACD": len(ACD), "ABD": len(ABD), "BCD": len(BCD),
                "ABCD": len(ABCD)}

            # don't bother correcting the scores yet, just do the venn:
            # a few corractions:
            # ABCD is subtracted from all other elements:

            # sub twos from singles:
            scores["A"] -= scores["AB"]
            scores["A"] -= scores["AC"]
            scores["B"] -= scores["AB"]
            scores["B"] -= scores["BD"]
            scores["C"] -= scores["CD"]
            scores["C"] -= scores["AC"]
            scores["D"] -= scores["CD"]
            scores["D"] -= scores["BD"]

            # then subtract the threes from their children:
            scores["AC"] -= scores["ABC"]
            scores["AC"] -= scores["ACD"]

            scores["AB"] -= scores["ABC"]
            scores["AB"] -= scores["ABD"]

            scores["CD"] -= scores["ACD"]
            scores["CD"] -= scores["BCD"]

            scores["BD"] -= scores["ABD"]
            scores["BD"] -= scores["BCD"]


            for k in scores:
                if k != "ABCD":
                    scores[k] -= scores["ABCD"]

            realfilename = self.draw._venn4(lists=lists, scores=scores, filename=filename)

        config.log.info("Saved Figure to '%s'" % realfilename)
        return(None)

    def moving_average_maps(self, mode="graph", compare_array=None, filename=None, key=None,
        normalise=True, **kargs):
        """
        **Purpose**
            Draw moving average maps in a variety of formats.

        **Arguments**
            mode (Optional, defaults to "graph")
                "graph"
                    draw a series of line graphs on the same graph.
                "heatmap"
                    draw as a heatmap
                "stacked_plots"
                    draw as a series of stacked line plots.

            compare_array (Required)
                The name of the array to use as a compare array. You have to draw the
                moving average by comparing against some array already in
                the glglob.

            filename
                filename to save the image as.

            key
                key to use to match between the compare_array and the rest of the
                data in this glglob

            normalise (Optional, default=True)
                normalise the data, True or False.

        **Returns**
            The actual filename used to save and an image saved
            to 'filename'
        """
        # get the compare array:
        compare_array = self[compare_array]
        assert compare_array, "the compare array was not found in this glglob"

        res = {}
        last_val = 0
        for pl in self.linearData:
            if pl.name != compare_array.name: # a simple == will fail here, so I use the names to compare
                res[pl.name] = []
                last_val = 0
                names = pl[key] # get the name column.
                for i, v in enumerate(compare_array):
                    if v[key] in names:
                        #last_val += 1
                        #res[pl.name].append(last_val)
                        res[pl.name].append(1)
                    else:
                        res[pl.name].append(0)

                res[pl.name] = utils.movingAverage(res[pl.name], int(len(res[pl.name]) * 0.2))[1] # keep only y,
                typical_length = len(res[pl.name]) # measure length to generate x axis later.

        # append the compare_array for comparison.

        res[compare_array.name] = [i[0] for i in compare_array["conditions"]][0:typical_length] # this will be a list of lists (with 1 entry) - flatten the list.
        # the arange fakes the x axis for plots.
        # I have to slice the array as it's slightly shorter than the movingAverage plots...
        # You will only notice this if the list is short and you can see the names

        if normalise: #this will normalise the y-axis
            for k in res:
                # normalise to 0-->100
                min_val = min(res[k])
                max_val = max(res[k]) - min_val

                for i, v in enumerate(res[k]):
                    res[k][i] = ((v-min_val) / max_val) * 100.0

        if mode == "graph":
            # fake the xdata:
            xdata = arange(0,typical_length)
            plot.cla()
            fig = plot.figure(figsize=(8,5))
            axis = fig.add_subplot(111)
            for pl in res:
                axis.plot(xdata, res[pl][1], label=pl)

            axis.set_title("")
            axis.legend(loc=2, markerscale=0.1)#, prop={"size": "xx-small"})
            fig.savefig(filename)
            real_filename = filename
        elif mode == "heatmap":
            # At the moment the data is in the form: [ [ (x,y), (x,y), ]:class ... []]
            # I need to strip out the x data.
            scalebar_name = "score"
            if normalise:
                scalebar_name = "normalised score"
            # use key to build an array ready for draw._heatmap
            colnames=[]
            for k in res:
                res[k] = res[k]
                colnames.append(k)

            real_filename = self.draw._heatmap(data=res, filename=filename, col_names=colnames, row_names=compare_array["name"],
                row_cluster=False, col_cluster=True, colour_map=cm.Blues, #vmax=1,
                colbar_label=scalebar_name)
        config.log.info("Saved image to '%s'" % real_filename)
