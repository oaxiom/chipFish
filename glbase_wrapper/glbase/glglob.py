"""

glglobs are almagamated lists, useful for drawing comparisons between lists.

glglobs also inherit from genelist, but many of the methods are broken.

renamed to glglob as it clashes with a matplotlib module name

"""

from __future__ import division

import sys, os, csv, string, math

from numpy import array, zeros
from array import array as qarray
from copy import deepcopy

import config
from flags import *
from genelist import genelist
from draw import draw
from errors import AssertionError

import matplotlib.pyplot as plot
import matplotlib.cm as cm
from scipy.stats import spearmanr, pearsonr

class glglob(genelist): # cannot be a genelist, as it has no keys...
    def __init__(self, *args, **kargs):
        genelist.__init__(self)

        # args should be a list of lists.
        # we then store them in the linearData set
        for a in args:
            for b in args:
                assert a == b, "Lists must be equivalent"

        self.linearData = args
        self._optimiseData()

    def __repr__(self):
        return("glbase.glglob")

    def __str__(self):
        return("glglob contains: %s items" % len(self.linearData))

    def _optimiseData(self): # no keys, so would die.
        pass

    def loadCSV(self):
        print "Error: glglobs cannot be sensibly represented as a CSV file, use load() to save binary copies"
        return(False)

    def saveCSV(self):
        print "Error: glglobs cannot be sensibly represented as a CSV file, use save() to save binary copies"
        return(False)

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
                No other methods are currently implemented.

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
        assert method in valid_methods, "must use a valid method for comparison (%s)" % ", ".join(valid_methods)
        assert filename, "Filename must be valid"
        assert key in self.linearData[0].linearData[0], "key '%s' not found" % (key) # just check one of the lists

        if not config.SILENT: print "Info: This may take a while, all lists are intersected by '%s' with '%s' key" % (method, key)

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

        if not config.SILENT: print "Info: Saved Figure to '%s'" % realfilename
        return(result_table)

if __name__ == "__main__":
        from peaklist import peaklist
        import numpy
        # get some data;
        data1 = peaklist(filename="tests/testA.csv")
        data2 = peaklist(filename="tests/testB.csv")
        data3 = peaklist(filename="tests/testC.csv")
        data4 = peaklist(filename="example/ccat_list.region", format=format_ccat_output)
        g = glglob(data1, data2, data3, data4, type="peaklist")

        r = g.compare(key="loc", filename="matrix_compare.png", method="collide", distance="euclidean")
