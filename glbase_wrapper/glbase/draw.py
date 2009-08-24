"""
draw class for glbase

this is a static class containing various generic methods for drawing etc.

It should be overridable with other draw classes and should be
drawing platform agnostic - at the moment it uses just matplotlib,
but in the future that may change. And when glbase is part of
chipFish it may use Cairo instead for a much more complex interactive
drawing session.

**TODO**

* move all of the drawing stuff into here.
* draw must be silent on all _ prefixed methods.
* all args should be **kargs.
* problems caused by embedding pylab into draw. It breaks copyability.
  And is probably not desirable anyways.
* unify _plot and _qplot

"""

import copy

from numpy import array, arange
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from scipy.spatial.distance import pdist # not in scipy.cluster.hierarchy.distance as you might expect :(
import numpy as np
from array import array as qarray

import config
from flags import *
from errors import AssertionError

if config.MATPLOTLIB_AVAIL:
    import matplotlib.pyplot as plot
    import matplotlib.cm as cm

# define static class here.
class draw:
    def __init__(self, parent):
        """
        draw class
        contains generic drawing procedures.
        please wrap your own around the _prefixed versions.
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

    def _heatmap(self, **kargs):
        """
        my own version of heatmap.

        This will draw a dendrogram ... etc...

        See the inplace variants as to how to use.
        """
        assert kargs.has_key("filename"), "_heatmap missing filename"
        assert kargs.has_key("data"), "_heatmap missing data"

        # default options and acceptable kargs:
        row_cluster = True
        col_cluster = True
        dpi = config.DEFAULT_DPI
        vmin = 0 # colour min max.
        vmax = 1

        # The data key should be a serialised Dict, I need to make an array.
        data = array([kargs["data"][key] for key in kargs["data"]])

        # deal with arguments;

        if kargs.has_key("bracket"):
            data = self._bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]
        if kargs.has_key("row_cluster"): row_cluster = kargs["row_cluster"]
        if kargs.has_key("dpi"): dpi = kargs["dpi"]
        if kargs.has_key("col_cluster"): col_cluster = kargs["col_cluster"]

        # a few grace and sanity checks here;
        if len(data) == 1: col_cluster = False # clustering with a single point?
        if len(data[0]) == 1: row_cluster = False # ditto.

        fig = plot.figure(dpi=dpi)

        if row_cluster:
            # ---------------- Fist plot (tree) ----------------------------
            ax1 = fig.add_subplot(141)

            # from scipy;
            # generate the dendrogram
            Y = pdist(data.T, metric='euclidean')
            Z = linkage(Y, 'single')
            a = dendrogram(Z, orientation='right')

            ax1.set_position([0.05,0.05,0.248,0.88])
            ax1.set_frame_on(False)
            ax1.set_xticklabels("")
            ax1.set_yticklabels("")
            ax1.set_ylabel("")
            # clear the ticks.
            tick_lines = ax1.yaxis.get_ticklines() # hack to blank the tick lines.
            for item in ax1.yaxis.get_ticklines():
                item.set_markeredgewidth(0.0)

            tick_lines = ax1.xaxis.get_ticklines() # hack to blank the tick lines.
            for item in ax1.xaxis.get_ticklines():
                item.set_markeredgewidth(0.0)

            # USe the tree to reorder the data.
            order = a["ivl"]

            # resort the data by order;
            newd = []
            new_row_names = []
            for index in order:
                newd.append(data.T[int(index)])
                new_row_names.append(kargs["row_names"][int(index)])
            data = array(newd)
            kargs["row_names"] = new_row_names
            # resort the row_names:
        else:
            data = data.T

        if col_cluster:
            ax2 = fig.add_subplot(142) # colour bar goes in here.
            ax2.set_frame_on(False)
            ax2.set_position([0.3,0.91,0.45,0.09])
            Y = pdist(data.T, metric='euclidean')
            Z = linkage(Y, 'single')
            a = dendrogram(Z, orientation='top')

            for item in ax2.yaxis.get_ticklines():
                item.set_markeredgewidth(0.0)
            for item in ax2.xaxis.get_ticklines():
                item.set_markeredgewidth(0.0)
            ax2.set_xticklabels("")
            ax2.set_yticklabels("")

            order = a["ivl"]
            # resort the data by order;
            newd = []
            new_col_names = []
            for index in order:
                newd.append(data.T[int(index)])
                new_col_names.append(kargs["col_names"][int(index)])
            data = array(newd).T
            kargs["col_names"] = new_col_names
        else:
            data = data

        kargs["data"] = data
        # ---------------- Second plot (heatmap) -----------------------
        ax3 = fig.add_subplot(143)
        hm = ax3.pcolor(data, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax)

        ax3.set_frame_on(False)
        ax3.set_position([0.3,0.05,0.45,0.88])
        ax3.set_xticks(arange(len(kargs["col_names"]))+0.5)
        ax3.set_xticklabels(kargs["col_names"])
        ax3.set_yticks(arange(len(kargs["row_names"]))+0.5)
        ax3.set_yticklabels(kargs["row_names"])
        ax3.yaxis.tick_right()
        for item in ax3.yaxis.get_ticklines():
            item.set_markeredgewidth(0.0)
        for item in ax3.xaxis.get_ticklines():
            item.set_markeredgewidth(0.0)
        for t in ax3.get_yticklabels(): # generally has to go last.
            t.set_fontsize(6)

        ax0 = fig.add_subplot(144)
        ax0.set_position([0.01, 0.96, 0.24, 0.03])
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0)
        cb.set_label("expression")

        for label in ax0.get_xticklabels():
            label.set_fontsize(5)

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
                            if kargs.has_key("match_key"):match_key = kargs["match_key"]
                            else:
                                print "Error: No suitable matching key between the microarray and the peaklist"
                                print "       valid keys are: refseq entrez loc tss_loc array_systematic_name"
                                print "       Both the microarray and the peaklist must have both keys"
                                sys.quit()

        # ----------------------modify defaults:
        if kargs.has_key("window"): moving_window = kargs["window"]
        if kargs.has_key("match_key"): match_key = kargs["match_key"]
        if kargs.has_key("dpi"): dpi = kargs["dpi"]
        if kargs.has_key("bracket"):
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]
        # bracketing is done below.

        # sanity checking

        # get arrange and sort the data.
        array_dict_data = kargs["array_data"].serialisedArrayDataDict
        arraydata = array([array_dict_data[key] for key in array_dict_data]) # needs to be sorted already.

        if kargs.has_key("bracket"):
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
        for item in ax0.yaxis.get_ticklines():
            item.set_markeredgewidth(0.0)

        ax1 = fig.add_subplot(142)

        hm = ax1.pcolor(arraydata.T, cmap=cm.RdBu_r)

        # set the display preferences for ax1 (the heatmap)
        ax1.set_frame_on(False)
        ax1.set_position([0.10,0.05,0.40,0.85])
        # y
        ax1.set_yticks(arange(len(kargs["row_names"]))+0.5)
        ax1.set_yticklabels(kargs["row_names"])
        ax1.yaxis.tick_left()
        for label in ax1.get_yticklabels():
            label.set_fontsize(4)

        #x
        ax1.set_xticks(arange(len(kargs["col_names"]))+0.5)
        ax1.set_xticklabels(kargs["col_names"])
        for item in ax1.xaxis.get_ticklines():
            item.set_markeredgewidth(0.0)
        for item in ax1.yaxis.get_ticklines():
            item.set_markeredgewidth(0.0)
        for label in ax1.get_xticklabels():
            label.set_fontsize(6)

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

        for item in ax2.yaxis.get_ticklines():
            item.set_markeredgewidth(0.0)

        for item in ax2.xaxis.get_ticklines():
            item.set_markeredgewidth(0.0)

        # bargraph -----------------------------------------------------

        ax3 = fig.add_subplot(144)
        #ax3.plot(peakdata, arange(len(peakdata))) # doesn't use the movingAverage generated x, scale it across the entire graph.
        #print len(peakdata), len(x)
        ax3.plot(peakdata, x) # doesn't use the movingAverage generated x, scale it across the entire graph.
        ax3.set_frame_on(False)
        ax3.set_position([0.6,0.05,0.3,0.85])
        ax3.set_yticklabels("")
        ax3.set_xticklabels("")
        for item in ax3.yaxis.get_ticklines():
            item.set_markeredgewidth(0.0)
        for item in ax3.xaxis.get_ticklines():
            item.set_markeredgewidth(0.0)

        m = utils.mean(peakdata)
        s = utils.std(peakdata)

        ax3.axvline(x=m, color='gray', linestyle="-", linewidth=0.5)
        ax3.axvline(x=(m+s), color='r', linestyle=":", linewidth=0.5)

        return(self._saveFigure(fig, kargs["filename"], dpi=dpi))

    def _plot(self, filename, data, **kargs):
        """
        Internal very very thin wrapper around matplotlib's plot
        """
        dpi = config.DEFAULT_DPI
        if kargs.has_key("dpi"): dpi = kargs["dpi"]

        fig = plot.figure(dpi=dpi)
        axis = fig.add_subplot(111)
        if kargs.has_key("x"):
            axis.plot(kargs["x"], data)
        else:
            axis.plot(data)

        if kargs.has_key("title"): axis.set_title(kargs["title"])
        if kargs.has_key("xlabel"): axis.set_xlabel(kargs["xlabel"])
        if kargs.has_key("ylabel"): axis.set_ylabel(kargs["ylabel"])

        return(self._saveFigure(fig, filename, dpi=dpi))

    def _qplot(self, list_of_tuples_data, filename=None, labels=None, **kargs):
        """
        thin wrapper around plot.
        expects list_of_tuples_data to be a list of tuples containing X and Y data.
        # valid kargs:
        labels = labels for each line (a list or iterable)
        dpi = dpi. Uses config.DEFAULT_DPI if not specified
        """
        assert filename, "Internal: _qplot no filename"
        assert labels, "Internal: _qplot, labels must be a list"

        dpi = config.DEFAULT_DPI
        if kargs.has_key("dpi"): dpi = kargs["dpi"]

        plot.cla() # paranoia!
        fig = plot.figure(dpi=dpi)
        axis = fig.add_subplot(111)
        for index, item in enumerate(list_of_tuples_data):
            if index > 0:
                label = "random"
            axis.plot(item[0], item[1], label=labels[index])
        #self.draw.py.yaxis([0,max(item[1])])
        axis.legend()
        return(self._saveFigure(fig, filename))

    def _saveFigure(self, fig, filename, dpi=config.DEFAULT_DPI):
        valid_output_modes = ["png", "ps", "eps", "svg"]
        
        assert config.DEFAULT_DRAWER in valid_output_modes, "Error: '%s' is not a supported drawing mode" % config.DEFAULT_DRAWER

        save_name = "%s.%s" % ("".join(filename.split(".")[:-1]), config.DEFAULT_DRAWER)

        fig.savefig(save_name)
        return(save_name)
