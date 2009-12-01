"""
draw class for glbase

this is a static class containing various generic methods for drawing etc.

It should be overridable with other draw classes and should be
drawing platform agnostic - at the moment it uses just matplotlib,
but in the future that may change.

**TODO**

* move all of the drawing stuff into here. [Note done, see microarray]
* unify _heatmap_and_plot with _heatmap (Can it be doe at all?)

"""

import sys, os, copy

from numpy import array, arange, mean, max, min
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from scipy.spatial.distance import pdist # not in scipy.cluster.hierarchy.distance as you might expect :(
import numpy as np
from array import array as qarray
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from matplotlib.patches import Ellipse, Circle

import config
from flags import *
from errors import AssertionError, DrawNoMatchingKeys_Special

# this is a work around for a poor implementation of
# scipy.cluster.hierarchy.
# I may need to implement this myself later.
# This can deal with ~23,000 x 14 at least.
# No idea on the upper limit.
sys.setrecursionlimit(5000) # 5x larger recursion.

# define static class here.
class draw:
    def __init__(self, parent):
        """
        draw class
        contains generic-ish drawing procedures.
        please wrap your own around the _prefixed versions.

        There are no public methods here.
        """
        self.parent = parent
        self.output_mode = config.DEFAULT_DRAWER

    def _bracket_data(self, data, min, max):
        """
        brackets the data between min and max (ie. bounds the data with no scaling)

        This should be a helper?
        """
        ran = max - min
        newd = copy.deepcopy(data)
        for x, row in enumerate(data):
            for y, value in enumerate(row):
                if value < min:
                    newd[x][y] = min
                elif value > max:
                    newd[x][y] = max
        return(newd)

    def _heatmap(self, cluster_mode = "euclidean", row_cluster = True, col_cluster = True, dpi = config.DEFAULT_DPI,
        vmin = 0, vmax = None, colour_map = cm.RdBu_r, **kargs):
        """
        my own version of heatmap.

        This will draw a dendrogram ... etc...

        See the inplace variants as to how to use.

        row_names is very important as it describes the order of the data.
        cluster_mode = pdist method. = ["euclidean"] ??????!
        """
        assert kargs["filename"], "_heatmap() missing filename"

        colour_map = cm.RdBu_r

        # preprocess data
        if isinstance(kargs["data"], dict):
            # The data key should be a serialised Dict, I need to make an array.
            data = array([kargs["data"][key] for key in kargs["col_names"]]).T
        else:
            # the default is a numpy like array object which can be passed right through.
            data = kargs["data"]

        # positions of the items in the plot:
        left_side_tree =    [0.05,  0.10,   0.248,  0.80]
        top_side_tree =     [0.3,   0.902,   0.45,   0.09]
        heatmap_location =  [0.3,   0.10,   0.45,   0.80]
        scalebar_location = [0.01,  0.96,   0.24,   0.03]
        row_font_size = 6
        col_font_size = 6

        if "square" in kargs and kargs["square"]:
            # make the heatmpa square, for e.g. comparison plots
            left_side_tree =    [0.198,  0.10,   0.10,   0.75]
            top_side_tree =     [0.3,   0.852,  0.55,   0.10]
            heatmap_location =  [0.3,   0.10,   0.55,   0.75]
            row_font_size = 10
            col_font_size = 10

        if "bracket" in kargs:
            data = self._bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]

        if not vmax:
            """
            I must guess the vmax value. I will do this by working out the
            mean then determining a symmetric colour distribution
            """
            me = mean(data)
            ma = abs(me - max(data))
            mi = abs(min(data) + me)
            if ma > mi:
                vmin = me - ma
                vmax = me + ma
            else:
                vmin = me - mi
                vmax = me + mi

        if "row_cluster" in kargs: row_cluster = kargs["row_cluster"]
        if "dpi" in kargs: dpi = kargs["dpi"]
        if "col_cluster" in kargs: col_cluster = kargs["col_cluster"]
        if not "colbar_label" in kargs: kargs["colbar_label"] = "expression"
        if "cmap" in kargs: colour_map = kargs["cmap"]

        # a few grace and sanity checks here;
        if len(data) <= 1: row_cluster = False # clustering with a single point?
        if len(data[0]) <= 1: col_cluster = False # ditto.

        fig = plot.figure(dpi=dpi)

        #print kargs["data"]
        #print data

        if row_cluster:
            # ---------------- Left side plot (tree) -------------------
            ax1 = fig.add_subplot(141)

            # from scipy;
            # generate the dendrogram
            Y = pdist(data, metric=cluster_mode)
            Z = linkage(Y, 'single')
            a = dendrogram(Z, orientation='right')

            ax1.set_position(left_side_tree)
            ax1.set_frame_on(False)
            ax1.set_xticklabels("")
            ax1.set_yticklabels("")
            ax1.set_ylabel("")
            # clear the ticks.
            [item.set_markeredgewidth(0.0) for item in ax1.xaxis.get_ticklines()]
            [item.set_markeredgewidth(0.0) for item in ax1.yaxis.get_ticklines()]

            # USe the tree to reorder the data.
            order = a["ivl"]

            # resort the data by order;
            if kargs["row_names"]: # make it possible to cluster without names
                newd = []
                new_row_names = []
                for index in order:
                    newd.append(data[int(index)])
                    new_row_names.append(kargs["row_names"][int(index)])
                data = array(newd)
                kargs["row_names"] = new_row_names
        else:
            data = data

        if col_cluster:
            # ---------------- top side plot (tree) --------------------
            transposed_data = data.T

            ax2 = fig.add_subplot(142)
            ax2.set_frame_on(False)
            ax2.set_position(top_side_tree)
            Y = pdist(transposed_data, metric=cluster_mode)
            Z = linkage(Y, 'single')
            a = dendrogram(Z, orientation='top')

            [item.set_markeredgewidth(0.0) for item in ax2.xaxis.get_ticklines()]
            [item.set_markeredgewidth(0.0) for item in ax2.yaxis.get_ticklines()]
            ax2.set_xticklabels("")
            ax2.set_yticklabels("")

            order = a["ivl"]
            # resort the data by order;
            if kargs["col_names"]: # make it possible to cluster without names
                newd = []
                new_col_names = []
                for index in order:
                    newd.append(transposed_data[int(index)])
                    new_col_names.append(kargs["col_names"][int(index)])
                data = array(newd).T # transpose back orientation
                kargs["col_names"] = new_col_names
        else:
            data = data

        # ---------------- Second plot (heatmap) -----------------------
        ax3 = fig.add_subplot(143)
        hm = ax3.pcolor(data, cmap=colour_map, vmin=vmin, vmax=vmax)

        ax3.set_frame_on(False)
        ax3.set_position(heatmap_location)
        if kargs["col_names"]:
            ax3.set_xticks(arange(len(kargs["col_names"]))+0.5)
            ax3.set_xticklabels(kargs["col_names"])
            ax3.set_xlim([0, len(kargs["col_names"])])
        else:
            ax3.set_xlim([0,data.shape[1]])

        if "square" in kargs and kargs["square"]:
            ax3.set_xticklabels(kargs["col_names"], rotation="vertical")

        if kargs["row_names"]:
            ax3.set_yticks(arange(len(kargs["row_names"]))+0.5)
            ax3.set_ylim([0, len(kargs["row_names"])])
            ax3.set_yticklabels(kargs["row_names"])
        else:
            ax3.set_ylim([0,data.shape[0]])

        ax3.yaxis.tick_right()
        [item.set_markeredgewidth(0.0) for item in ax3.xaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax3.yaxis.get_ticklines()]
        [t.set_fontsize(row_font_size) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(col_font_size) for t in ax3.get_xticklabels()]

        ax0 = fig.add_subplot(144)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label(kargs["colbar_label"])
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]

        return(self._saveFigure(fig, kargs["filename"], dpi=dpi))

    def _heatmap_and_plot(self, **kargs):
        """
        a two panel dendrogram - heatmap and plot.

        todo:

        * make it more generic - should be heatmap_data, plot_data.
        * some sort of helper ot show the interaction between the data sets.

        valid args:

        (Required)

        filename

            the filename to save the png to.

        array_data

            the data, usually serialisedArrayDataDict

        peak_data

            locations of the TF binding sites.

        (Optional)

        window

            the size of the moving average window (defaults to 50)

        match_key

            the key to match between the array and the peaklist.
            (Goes throgh a cascade until it finds one of:
            "refseq" > "entrez" > "loc" > "tss_loc" > "array_systematic_name")
        """

        # ----------------------defaults:
        moving_window = 50
        dpi = config.DEFAULT_DPI
        vmin = 0.0
        vmax = 1.0
        # work out a suitable match_key
        keys1 = kargs["array_data"].getKeys()
        keys2 = kargs["peak_data"].getKeys()
        # one day, I will think of a more elegant way to do this...
        # bitwise OR, then
        # this should be in peaklist, not here...
        if "refseq" in keys1 and "refseq" in keys2:
            match_key = "refseq"
        else:
            if "entrez" in keys1 and "entrez" in keys2:
                match_key = "entrez"
            else:
                if "loc" in keys1 and "loc" in keys2:
                    match_key = "loc"
                else:
                    if "tss_loc" in keys1 and "tss_loc" in keys2:
                        match_key = "tss_loc"
                    else:
                        if "array_systematic_name" in keys1 and "array_systematic_name" in keys2:
                            match_key = "array_systematic_name"
                        else:
                            if "match_key" in kargs:
                                match_key = kargs["match_key"]
                            else:
                                raise DrawNoMatchingKeys_Special

        # ----------------------modify defaults:
        if "window" in kargs: moving_window = kargs["window"]
        if "match_key" in kargs: match_key = kargs["match_key"]
        if "dpi" in kargs: dpi = kargs["dpi"]
        if "bracket" in kargs:
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]
        # bracketing is done below.

        # get arrange and sort the data.
        array_dict_data = kargs["array_data"].serialisedArrayDataDict
        arraydata = array([array_dict_data[key] for key in array_dict_data]) # needs to be sorted already.
        # There is a potential bug here, with the column names if multiple data is sent.

        if "bracket" in kargs:
            arraydata = self._bracket_data(arraydata, kargs["bracket"][0], kargs["bracket"][1])

        peakdata = kargs["peak_data"]

        bin = [0 for x in arraydata.T]

        for index, item in enumerate(kargs["array_data"]):
            f = peakdata._findByLabel(match_key, item[match_key])
            if f:
                bin[index] = 1

        x, peakdata = utils.movingAverage(bin, moving_window)

        # Now do the plots:
        plot.subplot(111)
        plot.cla()

        fig = plot.figure(dpi=dpi)

        # heatmap ------------------------------------------------------

        ax0 = fig.add_subplot(141) # colour bar goes in here.
        ax0.set_frame_on(False)
        ax0.set_position([0.10, 0.97, 0.30, 0.02])
        [item.set_markeredgewidth(0.0) for item in ax0.yaxis.get_ticklines()]


        ax1 = fig.add_subplot(142)

        hm = ax1.pcolor(arraydata.T, cmap=cm.RdBu_r)

        # set the display preferences for ax1 (the heatmap)
        ax1.set_frame_on(False)
        ax1.set_position([0.10,0.05,0.40,0.85])
        # y
        ax1.set_yticks(arange(len(kargs["row_names"]))+0.5)
        ax1.set_yticklabels(kargs["row_names"])
        ax1.yaxis.tick_left()
        [label.set_fontsize(4) for label in ax1.get_yticklabels()]

        #x
        ax1.set_xticks(arange(len(kargs["col_names"]))+0.5)
        ax1.set_xticklabels(kargs["col_names"])
        [item.set_markeredgewidth(0.0) for item in ax1.xaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax1.yaxis.get_ticklines()]
        [label.set_fontsize(6) for label in ax1.get_xticklabels()]

        hm = ax1.pcolor(arraydata.T, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)

        fig.colorbar(hm, cax=ax0, orientation="horizontal")
        for label in ax0.get_xticklabels():
            label.set_fontsize(6)

        # binding map --------------------------------------------------

        ax2 = fig.add_subplot(143)

        a = array(bin) # reshape the bin array
        a.shape = 1,len(bin)
        ax2.pcolor(a.T, cmap=cm.gray_r)

        ax2.set_frame_on(False)
        ax2.set_position([0.52,0.05,0.04,0.85])
        ax2.set_yticks(arange(len(kargs["row_names"]))+0.5)
        ax2.set_yticklabels("")
        ax2.set_xticklabels("")
        ax2.yaxis.tick_left()
        [item.set_markeredgewidth(0.0) for item in ax2.yaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax2.xaxis.get_ticklines()]

        # bargraph -----------------------------------------------------

        ax3 = fig.add_subplot(144)
        #ax3.plot(peakdata, arange(len(peakdata))) # doesn't use the movingAverage generated x, scale it across the entire graph.
        #print len(peakdata), len(x)
        ax3.plot(peakdata, x) # doesn't use the movingAverage generated x, scale it across the entire graph.
        ax3.set_frame_on(False)
        ax3.set_position([0.6,0.05,0.3,0.85])
        ax3.set_yticklabels("")
        ax3.set_xticklabels("")
        ax3.set_xlim([0, len(kargs["col_names"])])
        ax3.set_ylim([0, len(kargs["row_names"])])
        [item.set_markeredgewidth(0.0) for item in ax3.yaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax3.xaxis.get_ticklines()]

        m = utils.mean(peakdata)
        s = utils.std(peakdata)

        ax3.axvline(x=m, color='gray', linestyle="-", linewidth=0.5)
        ax3.axvline(x=(m+s), color='r', linestyle=":", linewidth=0.5)

        return(self._saveFigure(fig, kargs["filename"], dpi=dpi))

    def _boxplot(self, data=None, filename=None, **kargs):
        """
        super thin wrapper around matplotlib's boxplot
        """
        assert data, "data not found"
        assert filename, "no filename specified"

        dpi = config.DEFAULT_DPI
        if "dpi" in kargs: dpi = kargs["dpi"]

        plot.cla()
        fig = plot.figure(dpi=dpi)
        fig.clear()
        axis = fig.add_subplot(111)
        axis.boxplot(data)

        # handlers to modify common matplotlib options:
        if "xticklabels" in kargs:
            axis.set_xticklabels(kargs["xticklabels"])

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylable" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])

        return(self._saveFigure(fig, filename, dpi=dpi))

    def _scatter(self, x=None, y=None, filename=None, **kargs):
        """
        super thin wrapper aroung matplotlib's scatter
        """
        assert len(x), "x data missing"
        assert len(y), "y data missing"
        assert filename, "filename missing"

        dpi = config.DEFAULT_DPI
        if "dpi" in kargs:
            dpi = kargs["dpi"]

        plot.cla()
        fig = plot.figure(dpi=dpi)
        fig.clear()
        axis = fig.add_subplot(111)
        axis.scatter(x,y)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylable" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])

        return(self._saveFigure(fig, filename, dpi=dpi))

    def _qhist(self, filename=None, data=None, bins=60, **kargs):
        """
        Very thin wrapper around matplotlibs' hist
        """
        assert filename, "Internal Error: _qhist no filename"
        assert data, "Internal Error: _qhist missing data"

        dpi = config.DEFAULT_DPI
        if "dpi" in kargs:
            dpi = kargs["dpi"]

        plot.cla()
        fig = plot.figure(dpi=dpi)
        axis = fig.add_subplot(111)
        axis.hist(data, bins=bins, rwidth=0.8, facecolor='orange', alpha=0.7)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylable" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])

        return(self._saveFigure(fig, filename, dpi=dpi))

    def _plot(self, filename=None, data=None, **kargs):
        """
        Internal very very thin wrapper around matplotlib's plot
        """
        dpi = config.DEFAULT_DPI
        if "dpi" in kargs:
            dpi = kargs["dpi"]

        plot.cla()
        fig = plot.figure(dpi=dpi)
        axis = fig.add_subplot(111)
        if "x" in kargs:
            axis.plot(kargs["x"], data)
        else:
            axis.plot(data)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylable" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])

        return(self._saveFigure(fig, filename, dpi=dpi))

    def _qplotxy(self, list_of_tuples_data, filename=None, labels=None, **kargs):
        """
        thin wrapper around plot.
        expects list_of_tuples_data to be a list of tuples containing X and Y data:

        [ ([0..x], [0..y]), ([],[]) .. ([],[]) ] ???
        # valid kargs:
        labels = labels for each line (a list or iterable)
        dpi = dpi. Uses config.DEFAULT_DPI if not specified
        """
        assert filename, "Internal Error: _qplot no filename"
        assert labels, "Internal Error: _qplot, labels must be a list"

        dpi = config.DEFAULT_DPI
        if "dpi" in kargs:
            dpi = kargs["dpi"]

        # set up figure
        plot.cla() # paranoia!
        fig = plot.figure(dpi=dpi)
        axis = fig.add_subplot(111)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylable" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])

        for index, item in enumerate(list_of_tuples_data):
            if index > 0:
                label = "random"
            axis.plot(item[0], item[1], label=labels[index])
        #self.draw.py.yaxis([0,max(item[1])])
        if labels: axis.legend()
        return(self._saveFigure(fig, filename))

    def _vennDiagram2(self, left, right, overlap, filename=None,
        proportional=False, plot_sig=False, labels=None,
        **kargs):
        """
        draw a venn Diagram.

        **Arguments**

            (Required Arguments)

            left (relative)
                left - overlap. _vennDiagram2() does not correct for
                overlap.

            right (relative)
                right - overlap. _vennDiagram2() does not correct for
                overlap.

            overlap
                number of sites in common between left and right

            filename
                save name

            labels
                A dictionary in this form:
                {"left": "", "right": "", "title": ""}

            (Optional Arguments)

            proportional (True|False, default False)
                Use proportional circles suggesting the sizes of the overlaps

            plot_sig (True|False, default False)
                plot a binomial p-value test result for the overlap.

            simulations (100 by default)
                number of simulations to run

            world_size (defaults to sum of left, right, overlap)
                If your world_size is bigger you need to send it as the
                simulations will overestimate the overlaps.

            (Supported generic kargs)

            dpi
                dpi of the output file.

        **Result**
            saves a venn Diagram and returns the actual path to the saved file.

        **todo**
            support multiple overlaps.
            support for weighted circles.

        """
        # no valid_args specified here. You should know what you are doing...
        assert filename, "No filename specified!"
        assert labels, "No labels!"

        dpi = config.DEFAULT_DPI
        if "dpi" in kargs:
            dpi = kargs["dpi"]

        plot.cla()
        fig = plot.figure(dpi=dpi)
        axis = fig.add_subplot(111)
        axis.set_position([0.02, 0.02, 0.96, 0.96])

        artists = []
        if not proportional:
            artists.append(Circle((7, 10), 5.5, alpha=0.5, facecolor="#FFA200")) # (loc), size
            artists.append(Circle((13, 10), 5.5, alpha=0.5, facecolor="#104BA9"))

            artists.append(Circle((7, 10), 5.5, alpha=1.0, fill=False, lw=2)) # I draw the circle twice to give bold outer lines.
            artists.append(Circle((13, 10), 5.5, alpha=1.0, fill=False, lw=2))

            axis.text(5, 10, str(left), size=30, ha="center", va="center")
            axis.text(10, 10, str(overlap), size=25, ha="center", va="center")
            axis.text(15, 10, str(right), size=30, ha="center", va="center")
        else:
            max_size = 5.5
            max_attraction = max_size * 2
            min_size = 1.0
            total_score = float(left + right + overlap)
            left_w =  ((left + overlap) / total_score) * max_size
            left_w2 =  ((left) / total_score) * max_attraction
            right_w = ((right + overlap) / total_score) * max_size
            right_w2 = ((right)/ total_score) * max_attraction

            centre_w = (overlap / total_score)

            delta_dist = centre_w

            # sanity checking:
            #if delta_dist > max_attraction: delta_dist = max_attraction
            #if left_w < min_size: left_w = min_size
            #if right_w < min_size: right_w = min_size

            #print delta_dist

            artists.append(Circle((4.5+left_w2, 10), left_w, alpha=0.5, facecolor="#FFA200")) # (loc), size
            artists.append(Circle((15.5-right_w2, 10), right_w, alpha=0.5, facecolor="#104BA9"))

            #artists.append(Circle((6, 10), left_weight, alpha=1.0, fill=False, lw=2)) # I draw the circle twice to give bold outer lines.
            #artists.append(Circle((14, 10), right_weight, alpha=1.0, fill=False, lw=2))

        for a in artists: # add all artists...
            axis.add_artist(a)

        axis.set_xticks([0,20])
        axis.set_yticks([0,20])
        # clear frame and axis markers:
        axis.set_frame_on(False)
        axis.set_yticklabels("")
        axis.set_xticklabels("")
        [i.set_markeredgewidth(0.0) for i in axis.yaxis.get_ticklines()]
        [i.set_markeredgewidth(0.0) for i in axis.xaxis.get_ticklines()]

        # add the labels:
        axis.text(5, 17, labels["left"], size=15, ha="center", va="center")
        axis.text(15, 17, labels["right"], size=15, ha="center", va="center")
        axis.text(10, 19, labels["title"], size=28, ha="center", va="center")

        return(self._saveFigure(fig, filename))

    def _vennDiagram3():
        """
        same as _vennDiagram2 except it performs a triple overlap.
        """
        raise NotImplementedError

    # _vennDiagram4?

    def _saveFigure(self, fig, filename, dpi=config.DEFAULT_DPI):
        valid_output_modes = ["png", "ps", "eps", "svg"]

        assert config.DEFAULT_DRAWER in valid_output_modes, "'%s' is not a supported drawing mode" % config.DEFAULT_DRAWER

        save_name = "%s.%s" % ("".join(filename.split(".")[:-1]), config.DEFAULT_DRAWER)

        fig.savefig(save_name)
        return(save_name)
