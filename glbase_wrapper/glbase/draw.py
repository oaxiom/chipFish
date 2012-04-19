"""
draw class for glbase

this is a static class containing various generic methods for drawing etc.

**TODO**

There's a change in direction here. Instead of draw containing lots of generic draw functions
instead its more like a set of wrappers around common ways to do matplotlib stuff.
Drawing inside genelists is fine, as long as it follows this paradigm:

fig = self.draw.getfigure(**kargs)

ax = ...

etc ...

self.draw.do_common_args(fig, **kargs)
filename = fig.savefigure(fig, filename)

It would probably be an improvement if the class part was removed.

Instead a series of methods, exposed by draw.method() at the module level would be better.

This makes them more like helpers for matplotlib than a full fledged object.

This could be easily refactored by changing lines like::

class genelist:
    ... init
        self.draw = draw()
        
        to 
        
        self.draw = draw
        
For now, until I refactor the code to remove lines like that.
Also I want to rename this file gldraw to remove name clashes.

Then it can go::

    gldraw.heatmap()
    gldraw.scatter()

"""

import sys, os, copy

from numpy import array, arange, mean, max, min, std, float32
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from scipy.spatial.distance import pdist # not in scipy.cluster.hierarchy.distance as you might expect :(
from scipy import polyfit, polyval
from scipy.stats import linregress
import numpy as np
from array import array as qarray # this is now replaced by numpy?
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from matplotlib.colors import ColorConverter, rgb2hex
import matplotlib.mlab as mlab
from matplotlib.patches import Ellipse, Circle

import config
from flags import *
from errors import AssertionError

# this is a work around in the implementation of
# scipy.cluster.hierarchy. It does some heavy
# recursion and even with relatively small samples quickly eats the
# available stack.
# I may need to implement this myself later.
# This can deal with ~23,000 x 14 at least.
# No idea on the upper limit.
sys.setrecursionlimit(5000) # 5x larger recursion.

# define static class here.
class draw:
    def __init__(self, bad_arg=None, **kargs):
        """please deprecate me"""
        pass

    def bracket_data(self, data, min, max):
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
        Scheduled for deprecation
        """
        return(self.heatmap(cluster_mode, row_cluster, col_cluster, dpi,
        vmin, vmax, cm.RdBu_r, **kargs))

    def heatmap(self, filename=None, cluster_mode="euclidean", row_cluster=True, col_cluster=True, 
        vmin=0, vmax=None, colour_map=cm.RdBu_r, col_norm=False, row_norm=False,
        **kargs):
        """
        my own version of heatmap.

        This will draw a dendrogram ... etc...

        See the inplace variants as to how to use.

        row_names is very important as it describes the order of the data.
        cluster_mode = pdist method. = ["euclidean"] ??????!
        
        **Arguments**
            data (Required)
                the data to use. Should be a 2D array for the heatmap.
                
            filename (Required)
                The filename to save the heatmap to.
                
            col_norm (Optional, default=False)
                normalise each column of data between 0 .. max => 0.0 .. 1.0
                
            row_norm (Optional, default=False)
                similar to the defauly output of heatmap.2 in R, rows are normalised 0 .. 1
                
            row_tree (Optional, default=False)
                provide your own tree for drawing. Should be a Scipy tree. row_labels and the data
                will be rearranged based on the tree, so don't rearrnge the data yourself
                
            row_font_size (Optional, default=guess suitable size)
                the size of the row labels (in points). If set this will also override the hiding of
                labels if there are too many elements. 
                
            col_font_size (Optional, default=8)
                the size of the column labels (in points)
                
        **Returns**
            The actual filename used to save the image.
        """
        assert filename, "heatmap() - no specified filename"

        # preprocess data
        if isinstance(kargs["data"], dict):
            # The data key should be a serialised Dict, I need to make an array.
            data = array([kargs["data"][key] for key in kargs["col_names"]]).T 
            # If the lists are not square then this makes a numpy array of lists. 
            # Then it will fail below with a strange error.
            # Let's check to make sure its square:
            ls = [len(kargs["data"][key]) for key in kargs["col_names"]]

            if not all(x == ls[0] for x in ls):  
                raise Exception, "Heatmap data not Square"
        else:
            # the default is a numpy like array object which can be passed right through.
            data = array(kargs["data"], dtype=float32)
        
        if col_norm:
            for col in xrange(data.shape[1]):   
                data[:,col] /= float(data[:,col].max())
        
        if row_norm:
            for row in xrange(data.shape[0]):
                mi = min(data[row,:])
                ma = max(data[row,:])
                data[row,:] = (data[row,:]-mi) / (ma-mi)

        # positions of the items in the plot:
        left_side_tree =    [0.05,  0.1,   0.248,  0.85]
        top_side_tree =     [0.3,   0.952,  0.25,   0.044]
        heatmap_location =  [0.3,   0.1,   0.25,   0.85]
        scalebar_location = [0.01,  0.96,   0.24,   0.03]
        
        # set size of the row text depending upon the number of items:
        if "row_font_size" in kargs:
            row_font_size = kargs["row_font_size"]         
        else:
            if "row_names" in kargs and kargs["row_names"]:
                if len(kargs["row_names"]) <= 100:
                    row_font_size = 8
                elif len(kargs["row_names"]) <= 200:
                    row_font_size = 3
                elif len(kargs["row_names"]) <= 300:
                    row_font_size = 2
                else:
                    config.log.warning("heatmap has too many row labels to be visible. Suppressing row_labels")
                    kargs["row_names"] = None
                    row_font_size = 1
                
        col_font_size = 8
        if "col_font_size" in kargs:
            col_font_size = kargs["col_font_size"]

        if "square" in kargs and kargs["square"]:
            # make the heatmap square, for e.g. comparison plots
            left_side_tree =    [0.15,    0.15,   0.10,   0.75]
            top_side_tree =     [0.25,    0.90,  0.55,   0.08]
            heatmap_location =  [0.25,    0.15,   0.55,   0.75]
            row_font_size = 10
            col_font_size = 10

        if "bracket" in kargs: # done here so clustering is performed on bracketed data
            data = self.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
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
        if "col_cluster" in kargs: col_cluster = kargs["col_cluster"]
        if not "colbar_label" in kargs: 
            kargs["colbar_label"] = "expression"
        if "cmap" in kargs: colour_map = kargs["cmap"]

        # a few grace and sanity checks here;
        if len(data) <= 1: row_cluster = False # clustering with a single point?
        if len(data[0]) <= 1: col_cluster = False # ditto.

        if not "aspect" in kargs:
            kargs["aspect"] = "long"
        fig = self.getfigure(**kargs)

        if row_cluster:
            # ---------------- Left side plot (tree) -------------------
            ax1 = fig.add_subplot(141)

            # from scipy;
            # generate the dendrogram
            if "row_tree" in kargs:
                Z = kargs["row_tree"]
            else:
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

            # Use the tree to reorder the data.
            order = a["ivl"]

            # resort the data by order;
            if "row_names" in kargs and kargs["row_names"]: # make it possible to cluster without names
                newd = []
                new_row_names = []
                for index in order:
                    newd.append(data[int(index)])
                    new_row_names.append(kargs["row_names"][int(index)])
                data = array(newd)
                kargs["row_names"] = new_row_names
            else: # no row_names, I still want to cluster
                newd = []
                for index in order:
                    newd.append(data[int(index)])
                data = array(newd)

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

        # ---------------- Second plot (heatmap) -----------------------
        ax3 = fig.add_subplot(143)
        hm = ax3.pcolor(data, cmap=colour_map, vmin=vmin, vmax=vmax, antialiased=False)

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
        ax3.set_xticklabels(kargs["col_names"], rotation="vertical")

        if "row_names" in kargs and kargs["row_names"]: # you can't meaningfully see >200 labels. So suppress them:
            ax3.set_yticks(arange(len(kargs["row_names"]))+0.5)
            ax3.set_ylim([0, len(kargs["row_names"])])
            ax3.set_yticklabels(kargs["row_names"])
        else:
            ax3.set_ylim([0,data.shape[0]])
            ax3.set_yticklabels("")

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

        return(self.savefigure(fig, filename))

    def _heatmap_and_plot(self, peak_data=None, match_key=None, 
        arraydata=None, peakdata=None, bin=None, **kargs):
        """
        a two panel dendrogram - heatmap and plot.

        todo:

        * make it more generic - should be heatmap_data, plot_data.
        * some sort of helper ot show the interaction between the data sets.

        valid args:

        Required:

        filename

            the filename to save the png to.

        array_data

            the data, usually serialisedArrayDataDict

        peak_data

            locations of the TF binding sites.

        match_key

            the key to match between the array and the peaklist.


        Optional:

        window

            the size of the moving average window (defaults to 10% of the list)
            If you specify a number it is the number of elements in the list
            and not a percentage.

        use_tag_score (defaults to False)

            use the key set by "use_tag_score" to determine the plot intesity for
            the peakdata. By default if the array and the peaklist
            match then it simply uses the
        """
        # ----------------------defaults:
        vmin = 0.0
        vmax = 1.0

        # ----------------------modify defaults:
        if "window" in kargs: moving_window = kargs["window"]
        if "match_key" in kargs: match_key = kargs["match_key"]
        if "bracket" in kargs:
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]

        # Positions of the items in the figure:
        left_heatmap = [0.10,  0.05,  0.20,  0.85]
        scale_bar =    [0.10,  0.97,  0.30,  0.02]
        binding_map =  [0.32,  0.05,  0.08,  0.85]
        freq_plot =    [0.42,   0.05,  0.4,   0.85]

        # Now do the plots:
        plot.subplot(111)
        plot.cla()

        fig = self.getfigure(**kargs)

        # heatmap ------------------------------------------------------

        ax0 = fig.add_subplot(141) # colour bar goes in here.
        ax0.set_frame_on(False)
        ax0.set_position(scale_bar)
        [item.set_markeredgewidth(0.0) for item in ax0.yaxis.get_ticklines()]

        ax1 = fig.add_subplot(142)
        plot_data = arraydata.T
        hm = ax1.pcolor(plot_data, cmap=cm.RdBu_r, vmin=vmin, vmax=vmax, antialiased=False)

        ax1.set_frame_on(False)
        ax1.set_position(left_heatmap)
        if "col_names" in kargs and kargs["col_names"]:
            ax1.set_xticks(arange(len(kargs["col_names"]))+0.5)
            ax1.set_xticklabels(kargs["col_names"])
            ax1.set_xlim([0, len(kargs["col_names"])])
        else:
            ax1.set_xlim([0,plot_data.shape[1]])
            ax1.set_xticklabels("")

        if "row_names" in kargs and kargs["row_names"] and len(kargs["row_names"]) <= 200: # you can't meaningfully see >200 labels. So suppress them:
            ax1.set_yticks(arange(len(kargs["row_names"]))+0.5)
            ax1.set_ylim([0, len(kargs["row_names"])])
            ax1.set_yticklabels(kargs["row_names"])
        else:
            ax1.set_ylim([0,plot_data.shape[0]])
            ax1.set_yticklabels("")

        ax1.yaxis.tick_left()
        [item.set_markeredgewidth(0.0) for item in ax1.xaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax1.yaxis.get_ticklines()]
        [t.set_fontsize(6) for t in ax1.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(6) for t in ax1.get_xticklabels()]

        fig.colorbar(hm, cax=ax0, orientation="horizontal")
        for label in ax0.get_xticklabels():
            label.set_fontsize(6)

        # binding map --------------------------------------------------

        ax2 = fig.add_subplot(143)

        a = array(bin) # reshape the bin array
        a.shape = 1,len(bin)
        ax2.pcolor(a.T, cmap=cm.gray_r, antialiased=False)

        ax2.set_frame_on(False)
        ax2.set_position(binding_map)
        ax2.set_yticks(arange(len(kargs["row_names"]))+0.5)
        ax2.set_yticklabels("")
        ax2.set_xticklabels("")
        ax2.set_xlim([0,1])
        ax2.set_ylim([0,len(kargs["row_names"])])
        ax2.yaxis.tick_left()
        [item.set_markeredgewidth(0.0) for item in ax2.yaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax2.xaxis.get_ticklines()]

        # linegraph -----------------------------------------------------

        ax3 = fig.add_subplot(144)
        ax3.plot(peakdata, arange(len(peakdata))) # doesn't use the movingAverage generated x, scale it across the entire graph.
        ax3.set_frame_on(False)
        ax3.set_position(freq_plot)
        ax3.set_yticklabels("")
        ax3.set_ylim([0, len(peakdata)])
        ax3.set_xlim([0, max(peakdata)])
        [item.set_markeredgewidth(0.0) for item in ax3.yaxis.get_ticklines()]
        [item.set_markeredgewidth(0.2) for item in ax3.xaxis.get_ticklines()]
        [t.set_fontsize(6) for t in ax3.get_xticklabels()]

        m = utils.mean(peakdata)
        s = utils.std(peakdata)

        ax3.axvline(x=m, color='black', linestyle=":", linewidth=1)
        ax3.axvline(x=(m+s), color='r', linestyle=":", linewidth=0.5)
        ax3.axvline(x=(m-s), color='r', linestyle=":", linewidth=0.5)

        return(self.savefigure(fig, kargs["filename"]))

    def boxplot(self, data=None, filename=None, labels=None, **kargs):
        """
        wrapper around matplotlib's boxplot
        """
        assert data, "data not found"
        assert filename, "no filename specified"
        assert labels, "boxplots must have labels"

        fig = self.getfigure(**kargs)

        ax = fig.add_subplot(111)
        r = ax.boxplot(data)

        plot.setp(r['medians'], color='red')
        plot.setp(r['whiskers'], color='black', lw=2)
        plot.setp(r['boxes'], color='black', lw=2)
        plot.setp(r['fliers'], color="grey")

        ax.set_xticklabels(labels)

        self.do_common_args(ax, **kargs)

        return(self.savefigure(fig, filename))

    def _scatter(self, x=None, y=None, filename=None, **kargs):
        """
        super thin wrapper aroung matplotlib's scatter
        """
        assert len(x), "x data missing"
        assert len(y), "y data missing"
        assert filename, "filename missing"

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)
        axis.scatter(x,y)

        if "logx" in kargs and kargs["logx"]: axis.set_xscale("log")
        if "logy" in kargs and kargs["logy"]: axis.set_yscale("log")

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylable" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])

        return(self.savefigure(fig, filename))

    def _qhist(self, filename=None, data=None, bins=60, **kargs):
        """
        Very thin wrapper around matplotlibs' hist
        """
        assert filename, "Internal Error: _qhist no filename"
        assert data, "Internal Error: _qhist missing data"

        log = False
        if "log" in kargs:
            log = kargs["log"]

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)
        axis.hist(data, bins=bins, facecolor='orange', ec="none", alpha=0.7, log=log)

        if "title" in kargs: 
            axis.set_title(kargs["title"])
        if "xlabel" in kargs: 
            axis.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs: 
            axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: 
            axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: 
            axis.set_yticks(kargs["yaxis"])

        return(self.savefigure(fig, filename))

    def _plot(self, filename=None, data=None, **kargs):
        """
        Internal very very thin wrapper around matplotlib's plot
        """

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)

        if "x" in kargs:
            axis.plot(kargs["x"], data)
        else:
            axis.plot(data)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])
        if "xlims" in kargs: axis.set_xlim(kargs["xlims"])
        if "ylims" in kargs: axis.set_ylim(kargs["ylims"])

        return(self.savefigure(fig, filename))

    def _genome_segment(self, figure_axis, loc, feature_list):
        """
        draw a representation of the genome using the axis provided by figure_axis
        loc = the genomic location
        feature_list = a list of refseq features to draw on the graph.
        """
        ax = figure_axis
        # set the x axis to match.
        ax.set_xlim([0, len(loc)]) # 1 = 1bp scale.
        ax.set_ylim([0, 10]) # arbitary.

        ax.set_xticks([0, len(loc)])
        ax.set_yticks([0, 10])

        [i.set_markeredgewidth(0.0) for i in ax.yaxis.get_ticklines()] # clear ticks
        [i.set_markeredgewidth(0.0) for i in ax.xaxis.get_ticklines()]
        left_base = loc["left"]
        for item in feature_list:
            if item["type"] == "gene":
                left_most_base = item["loc"]["left"] - loc["left"]
                right_most_base = item["loc"]["right"] - loc["left"]

                if item["strand"] in positive_strand_labels:
                    ax.text(left_most_base, 7, item["name"], size=10, ha="left", va="center")
                    ax.arrow(left_most_base, 6, right_most_base-left_most_base, 0,
                        alpha=0.5, fc=(0,0,0), width=0.02)
                elif item["strand"] in negative_strand_labels:
                    ax.text(right_most_base, 3, item["name"], size=10, ha="right", va="center")
                    ax.arrow(left_most_base, 4, right_most_base-left_most_base, 0,
                        alpha=0.5, fc=(0,0,0), width=0.02)
        ax.axhline(y=5, color='gray', linestyle=":", linewidth=0.5, alpha=0.5)
        # tidy up stuff:
        ax.set_frame_on(False)
        ax.set_yticklabels("")
        ax.set_xticklabels("")
        # append the chromosome coords to the figure.
        ax.text(100, 8, str(loc), size=10, ha="left", va="center")

    def _labeled_figure(self, data=None, axissize=None, filename=None,
        horizontal_line=True, **kargs):
        """
        **Purpose**
            Draw a figure with a set of labels.

        **Arguments**
            filename
                filename to save to. May get modified depending upon the current
                draw mode.

            data
                A set of values if this form: {"pos": (x,y), "label": label}

            horizontal_line (Optional, default=True)
                draw  horizontal line at y axis 0.

            axissize (Required)
                the axis dimensions (x and y maximal limits).

            Common Arguments also supported:
                figsize - tuple specifying the figure aspect
                dpi - the dpi (only supported for ps outputs)
        """
        position_plot = [0.02, 0.05, 0.96, 0.9]

        fig = self.getfigure(**kargs)

        ax1 = fig.add_subplot(131)
        ax1.set_position(position_plot)
        ax1.set_xlim([0, axissize[0]])
        ax1.set_ylim([0, axissize[1]])
        if "ylim" and kargs["ylim"]:
            ax1.set_ylim(kargs["ylim"])
        if horizontal_line:
            ax1.axhline(y=4.5, color='gray', linestyle=":", linewidth=1.5)

        for l in data:
            # draw an arrow.
            ax1.text(l["pos"][0], l["pos"][1], l["label"], size=8, ha="left", va="bottom",
                rotation="vertical")

        if "genomic_features" in kargs:
            # I've got some genomic features, I want to draw them.
            position_genomic = [0.02, 0.05, 0.96, 0.1]
            ax3 = fig.add_subplot(133)
            ax3.set_position(position_genomic)
            self._genome_segment(ax3, kargs["loc"], kargs["genomic_features"])

        return(self.savefigure(fig, filename))

    def _plot_and_histogram(self, filename=None, data=None, figsize=(5,5), **kargs):
        """
        Draw a graph plot and a histogram on the right hand side.
        """
        position_plot = [0.02, 0.05, 0.83, 0.9]
        position_histogram = [0.87, 0.05, 0.12, 0.9]


        fig = self.getfigure(**kargs)
        ax1 = fig.add_subplot(131)
        
        ax1.set_position(position_plot)
        if "x" in kargs: # fake some x data.
            ax1.plot(kargs["x"], data)
        else:
            ax1.plot(data)

        # histogram plot.
        ax2 = fig.add_subplot(132)
        ax2.set_position(position_histogram)
        n, bins, patches = ax2.hist(data, bins=20, orientation='horizontal', histtype="stepfilled", color=(0,0,0))

        m = mean(data) # hehe, numpy pawns my homemade routine for nontrivial samples.
        s = std(data)

        y = mlab.normpdf( bins, m, s)
        l = ax2.plot(bins, y, 'r--', linewidth=1)
        # add mean, std lines to the first plot.
        # add mean line
        ax1.axvline(x=m, color='black', linestyle="-", linewidth=0.5)
        # add std lines
        for z in [2,3,4]:
            zz = s*z # num stds away
            ax1.axhline(y=(m+zz), color='gray', linestyle=":", linewidth=1.5)
            ax1.axhline(y=(m-zz), color='gray', linestyle=":", linewidth=1)

        self.do_common_args(ax1, **kargs)
        if "xlims" in kargs: ax1.set_xlim(kargs["xlims"])
        if "ylims" in kargs: ax1.set_ylim(kargs["ylims"])

        #if "xlims" in kargs: ax2.set_xlim(kargs["xlims"])
        if "ylims" in kargs: ax2.set_ylim(kargs["ylims"])

        if "genomic_features" in kargs:
            # I've got some genomic features, I want to draw them.
            # I have to generate my own axis for draw._genome_segment
            position_genomic = [0.02, 0.05, 0.83, 0.1]
            ax3 = fig.add_subplot(133)
            ax3.set_position(position_genomic)
            self._genome_segment(ax3, kargs["loc"], kargs["genomic_features"])

        return(self.savefigure(fig, filename))

    # Salted for deprecation
    def _qplotxy(self, list_of_tuples_data, filename=None, labels=None, **kargs):
        self.qplotxy(list_of_tuples_data, filename=filename, labels=labels, **kargs)

    def qplotxy(self, list_of_tuples_data, filename=None, labels=None, **kargs):
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

        # set up figure
        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])
        if "xlim" in kargs: axis.set_xlim(kargs["xlim"])
        if "ylim" in kargs: axis.set_ylim(kargs["ylim"])

        for index, item in enumerate(list_of_tuples_data):
            axis.plot(item[0], item[1], label=labels[index])
            
        #self.draw.py.yaxis([0,max(item[1])])
        if labels: 
            leg = axis.legend()#bbox_to_anchor=(0., 1.02, 1., .102), loc=3)
            [t.set_fontsize(6) for t in leg.get_texts()]
            
        return(self.savefigure(fig, filename))

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

        fig = self.getfigure(**kargs)
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

        return(self.savefigure(fig, filename))

    def _venn3():
        """
        same as _vennDiagram2 except it performs a triple overlap.
        """
        raise NotImplementedError

    def _venn4(self, filename=None, lists=None, scores=None,
        proportional=False, **kargs):
        """
        draw a 4-way venn Diagram.

        **Arguments**
            matrix_grid (Required)
                A dictionary of the form:
                    {"A": A, "B": B, "C": D,
                    "AC": AB, "AB": AB, "CD": CD, "BD": BD,
                    "ABC": ABC, "ACD": ACD, "ABD": ABD, "BCD": BCD,
                    "ABCD": ABCD}
                containing the overlapping genelists
                (ie. the lists after a call to map()).

            The numbers should be corrected overlaps, (i.e,. A= A - AC - ABC - ABCD)

            filename
                save name

            (Optional Arguments)

            proportional (True|False, default False)
                EXPERIMENTAL!
                Use proportional circles suggesting the sizes of the overlaps

            NotImplemented:

            plot_sig (True|False, default False)
                plot a binomial p-value test result for the overlap.

            simulations (100 by default)
                number of simulations to run for the p-value

            world_size (defaults to sum of enitre data)
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


        m = lists
        s = scores

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)
        axis.set_position([0.02, 0.02, 0.96, 0.96])

        artists = []
        if not proportional:
            artists.append(Circle((8, 8), 5, alpha=0.7, facecolor="#9930a7")) # (loc), size


            artists.append(Circle((12, 8), 5, alpha=0.6, facecolor="#fd8348"))


            artists.append(Circle((12, 12), 5, alpha=0.5, facecolor="#30aa7f")) # (loc), size


            artists.append(Circle((8, 12), 5, alpha=0.4, facecolor="#e7f847"))

            # outlines:
            artists.append(Circle((8, 8), 5, alpha=1.0, fill=False, lw=2))
            artists.append(Circle((12, 8), 5, alpha=1.0, fill=False, lw=2))
            artists.append(Circle((12, 12), 5, alpha=1.0, fill=False, lw=2))
            artists.append(Circle((8, 12), 5, alpha=1.0, fill=False, lw=2)) # I draw the circle twice to give bold outer lines.


            # labels:
            axis.text(6, 6, "%s" %str(s["A"]), size=25, ha="center", va="center")
            axis.text(14, 6, "%s" %str(s["B"]), size=25, ha="center", va="center")
            axis.text(6, 14, "%s" %str(s["C"]), size=25, ha="center", va="center")
            axis.text(14, 14, "%s" %str(s["D"]), size=25, ha="center", va="center")

            # doubles:
            axis.text(10, 6, "%s" % str(s["AB"]), size=20, ha="center", va="center")
            axis.text(6,  10, "%s" % str(s["AC"]), size=20, ha="center", va="center")
            axis.text(10, 14, "%s" % str(s["CD"]), size=20, ha="center", va="center")
            axis.text(14, 10, "%s" % str(s["BD"]), size=20, ha="center", va="center")

            axis.text(12, 12, "%s" % str(s["BCD"]), size=15, ha="center", va="center")
            axis.text(8, 12, "%s" % str(s["ACD"]), size=15, ha="center", va="center")
            axis.text(8, 8,  "%s" % str(s["ABC"]), size=15, ha="center", va="center")
            axis.text(12, 8, "%s" % str(s["ABD"]), size=15, ha="center", va="center")

            # Fours:
            axis.text(10, 10, str(s["ABCD"]), size=25, ha="center", va="center")

            #axis.text(10, 10, str(overlap), size=25, ha="center", va="center")
            #axis.text(15, 10, str(right), size=30, ha="center", va="center")
        else:
            raise NotImplementedError

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
        axis.text(5, 2, m["A"].name, size=15, ha="center", va="center")
        axis.text(15, 2, m["B"].name, size=15, ha="center", va="center")
        axis.text(5, 18, m["C"].name, size=15, ha="center", va="center")
        axis.text(15, 18, m["D"].name, size=15, ha="center", va="center")
        #axis.text(15, 17, labels["right"], size=15, ha="center", va="center")
        #axis.text(10, 19, labels["title"], size=28, ha="center", va="center")

        return(self.savefigure(fig, filename))

    def getfigure(self, size=None, aspect=None, **kargs):
        """
        **Purpose**
            setup a valid figure instance based on size.
        
        **Arguments**
            size (Optional, default="medium")
                if size is a tuple then that tuple is the specified size in inches (don't ask)
                You can also specify "small", "medium", "large" and "huge". Corrsponding to approximate pixel
                sizes of (with the "normal" aspect)
                small  : 
                medium : 
                large  :
                huge   :
                
                If size is specified then it takes preference over aspect and aspect will be ignored.
                
            aspect (Optional, default="normal")
                the aspect of the image.
                currently only "normal", "long" and "square" are respected).
        
        **Returns**
            A valid matplotlib figure object
        """
        # options in the args seem to be set at compile time. 
        # So I have to interpret them here
        if not size:
            size = config.draw_size
        elif len(size) == 2: # A tuple or list?
            size_in_in = (size[0], size[1])

        if not aspect:
            aspect = config.draw_aspect
        
        if len(size) == 2:
            size_in_in = (size[0], size[1])
            return(plot.figure(figsize=size_in_in))
        else:
            data = {"normal": {"small": (5,4), "medium": (8,6), "large": (12,9), "huge": (16,12)},
                    "square": {"small": (4,4), "medium": (7,7), "large": (9,9), "huge": (12,12)},
                    "long": {"small": (4,5), "medium": (6,8), "large": (9,12), "huge": (12,16)},
                    "wide": {"small": (7,4), "medium": (12,6), "large": (18,9), "huge": (24,12)}
                    }
            dpi = {"small": 75, "medium": 150, "large": 300, "huge": 600} # This dpi doesn't actually work here...
            # See savefigure() for the actual specification
            return(plot.figure(figsize=data[aspect][size]))

    def savefigure(self, fig, filename, size=config.draw_size):
        """
        **Purpose**
            Save the figure
            to filename, modifying the filename based on the current drawing mode
            (if required)
        **Arguments**
            fig
                the figure handle
                
            filename
                the filename to save the file to
                
        **Returns**
            the actual filename used to save the image
        """
        assert config.draw_mode in config.valid_draw_modes, "'%s' is not a supported drawing mode" % config.draw_mode
        
        # So that saving supports relative paths.
        path, head = os.path.split(filename)
        save_name = "%s.%s" % (".".join(head.split(".")[:-1]), config.draw_mode) # this will delete .. in filename, e.g. file.meh.png
        
        dpi = {"small": 75, "medium": 150, "large": 200, "huge": 300}
        fig.savefig(os.path.join(path, save_name), dpi=dpi[size])
        plot.close(fig) # Saves a huge amount of memory.
        return(save_name)            

    def do_common_args(self, ax, **kargs):
        """
        **Purpose**
            deal with common arguments to matplotlib (may not always work, depending upon the figure type.
        
        **Arguments**
            ax
                an matplotlib axes object
        
            These are based loosly on the matplotlib versions
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits
                ylims - y-axis limits
                logx - set the x scale to a log scale argument should equal the base
                logy - set the y scale to a log scale
                legend_size - size of the legend, small, normal, medium
            
        **Returns**
            None
        """
        if "xlabel" in kargs:
            ax.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs:
            ax.set_ylabel(kargs["ylabel"])
        if "title" in kargs:
            ax.set_title(kargs["title"])
        if "xlims" in kargs:
            ax.set_xlim(kargs["xlims"])
        if "ylims" in kargs:
            ax.set_ylim(kargs["ylims"])      
        if "logx" in kargs:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs:
            ax.set_yscale("log", basey=kargs["logy"])
        if "legend_size" in kargs:
            legend = ax.get_legend()
            [t.set_fontsize(kargs["legend_size"]) for t in legend.get_texts()]
            

    def _stacked_plots(self, filename, loc, features, graphs, merged=False, **kargs):
        """
        Draw a stacked set of line plots

        graph_data is 'stackable', the more graphs added the graphs
        will seperate into different graphs.

        See pwm.scan_seq_with_features() for details of usage.

        graphs should be a dictionary with of which the key will be used as
        a name
        """

        dpi = config.DEFAULT_DPI
        if "dpi" in kargs:
            dpi = kargs["dpi"]

        fig = self.getfigure(**kargs)

        spacer = 0.9 / len(graphs)

        for index, key in enumerate(graphs):
            ax2 = fig.add_subplot(1, len(graphs), index)
            ax2.set_position([0.02, spacer*index, 0.96, spacer])
            ax2.plot(graphs[key])

            #ax2.set_frame_on(False)
            ax2.annotate(key, xy=(0.02, spacer*index))
            ax2.set_xticklabels("")
            ax2.set_yticklabels("")
            ax2.set_ylabel("")
            [item.set_markeredgewidth(0.0) for item in ax2.xaxis.get_ticklines()]
            [item.set_markeredgewidth(0.0) for item in ax2.yaxis.get_ticklines()]
            ax2.set_ylim([0, 1])
            ax2.set_xlim([0, len(graphs[key])])

        return(self.savefigure(fig, filename))
        
    def _simple_heatmap(self, filename=None, colour_map = cm.Reds, vmin=0, vmax=None, symmetric=False, **kargs):
        """
        A simplified version of heatmap, with no clustering, and a simpler representation
        Also, you can change the size and aspect of the display.
        
        **Arguments**
            data
                an array of arrays or equivalent.
        
            colour_map
                Default is YlOrRd
                
            bracket
                specify a tuple for the min max values for the heatmap.
                
            symmetric (Optional, default=False)
                If set to true, find the mean, the max and the min n the data
                and then set the range of colours to span min .. mean .. max, so that
                the gap between mean and max and min and max are identical.
                If set to False (default behaviour) then simply range the colours from min(data) to
                max(data).
                
            fig_size (Optional, default=(6,6))
                change the figure size aspect.
        """
        # This should be a wrapper around draw.heatmap() to take advantage of heatmaps 
        # better code.
        
        assert filename, "_heatmap() missing filename"

        data = kargs["data"]

        # positions of the items in the plot:
        heatmap_location =  [0.12,   0.01,   0.75,   0.98]
        scalebar_location = [0.01,  0.96,   0.10,   0.03]

        if "bracket" in kargs:
            data = self.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]

        if not vmax:
            """
            I must guess the vmax value. I will do this by working out the
            mean then determining a symmetric colour distribution
            """
            if symmetric:
                me = mean(data)
                ma = abs(me - max(data))
                mi = abs(min(data) + me)
                if ma > mi:
                    vmin = me - ma
                    vmax = me + ma
                else:
                    vmin = me - mi
                    vmax = me + mi
            else:
                vmax = max(data)
                vmin = min(data)
                
        if not "aspect" in kargs:
            kargs["aspect"] = "long"
        fig = self.getfigure(**kargs)

        # ---------------- Second plot (heatmap) -----------------------
        ax3 = fig.add_subplot(121)
        hm = ax3.pcolor(data, cmap=colour_map, vmin=vmin, vmax=vmax, antialiased=False)

        ax3.set_frame_on(True)
        ax3.set_position(heatmap_location)
        ax3.set_xlim([0,data.shape[1]])
        ax3.set_ylim([0,data.shape[0]])
        ax3.set_yticklabels("")
        ax3.set_xticklabels("")
        ax3.yaxis.tick_right()
        [item.set_markeredgewidth(0.0) for item in ax3.xaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax3.yaxis.get_ticklines()]
        #[t.set_fontsize(6) for t in ax3.get_yticklabels()] # generally has to go last.
        #[t.set_fontsize(6) for t in ax3.get_xticklabels()]

        ax0 = fig.add_subplot(122)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label("")
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]

        return(self.savefigure(fig, filename))
    
    def nice_scatter(self, x=None, y=None, filename=None, do_best_fit_line=False, 
        print_correlation=False, **kargs):
        """
        **Purpose**
            Draw a nice simple scatter plot
            
        **Arguments**
            x, y (Required)
                x and y values
                
            filename (Required)
                the filename to save as.
                
            spots (Optional, must be a 2-length tuple containing (x, y) data)
                These spots will be empahsised with whatever spots_cols is or an 
                "Orange" colour by default
                
            do_best_fit_line (Optional, default=False)
                Draw a line of best fit and the
                
            print_correlation (Optional, default=None)
                You have to spectify the type of correlation to print on the graph.
                valid are:
                    pearson = pearson correlation.
                    pearsonr2 = R^2.
                
                NOTE: do_best_fit_line must be True for this to work.
            
            Supported keyword arguments:
                xlabel, ylabel, title, logx, logy
        
        **Returns**
            the real filename, which may get modified depending upon the current drawing mode
            (usually results in a png)
        """
        fig = self.getfigure(aspect="square")
        ax = fig.add_subplot(111)
        
        ax.scatter(x, y, s=2, c="black", alpha=0.2, edgecolors="none")
        
        if "spots" in kargs and kargs["spots"]:
            if "spots_cols" in kargs and kargs["spots_cols"]:
                # Will recognise a string or sequence autmagivally.
                ax.scatter(kargs["spots"][0], kargs["spots"][1], s=5, c=kargs["spots_cols"], alpha=0.8, edgecolor="none")
            else:
                ax.scatter(kargs["spots"][0], kargs["spots"][1], s=5, c="orange", alpha=0.8, edgecolor="none")
        
        if do_best_fit_line:
            # linear regression
            (ar, br) = polyfit(x, y, 1)
            xr = polyval([ar,br], x)
            slope, intercept, r_value, p_value, std_err = linregress(x,y)
    
            mx = [min(x), max(x)]
            my = [slope * min(x) + intercept, slope * max(x) + intercept]
    
            ax.plot(mx, my, "r.-")
            
            if print_correlation:
                if print_correlation == "pearson":
                    ax.text(mx[0], my[0], "R=%.2f" % r_value)
                elif print_correlation == "pearsonr2":
                    ax.text(mx[0], my[0], "R2=%.2f" % (r_value*r_value))
        
        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])
        
        self.do_common_args(ax, **kargs)
        
        return(self.savefigure(fig, filename))
    
    def bar_chart(self, filename=None, genelist=None, data=None, cols=None, **kargs):
        """
        **Purpose**
            draw a bar chart with error bars and interpret and package the data coming from a genelist-like object
            
        **Args**
            filename
            
            genelist
                a genelist-like object
                
            data
                the key to look for in the genelist for the data
                
            labels
                the key to look for in the genelist for labels
                
            title (Optional)
                the title
                
            err (Optional)
                the key to look for in the genelist for error bar values.
                This one assumes symmetric values +- around the data
            
            err_up (Optional)
                the key to look for errorbars going up
                
            err_dn (Optional)
                the key to look for error bars going down
                
            errs_are_absolute (Optional, default=False)
                error bars are not +- from the data, but are values that specify where the error
                bars extend to. This needs to be set to True commonly for confidence intervals
                and left as False for standard errors.
            
            cols (Optional, default=Use a default set from matplotlib)
                the colours to use for the bar charts, there shoudl be one for each bar.
                
            Other kargs respected by bar_chart:
                aspect
                size
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits
                ylims - y-axis limits
                logx - set the x scale to a log scale argument should equal the base
                logy - set the y scale to a log scale
            
        **Returns**
            The actual_filename used to save the image
        """
        
        da = []
        err = []
        err_up = []
        err_dn = []
        for i in genelist:
            da.append(i[data])
            if "err" in kargs:
                err.append(i[kargs["errs"]])
            if "err_up" in kargs:
                err_up.append(i[kargs["err_up"]])
            if "err_dn" in kargs:
                err_dn.append(i[kargs["err_dn"]])
          
        if "errs_are_absolute" in kargs and kargs["errs_are_absolute"]:
            for i, n in enumerate(da): # normalise the values so matplotlib can underastand them
                if err_up:
                    err_up[i] = [a - b for a, b in zip(err_up[i], n)]
                if err_dn:
                    err_dn[i] = [b - a for a, b in zip(err_dn[i], n)]
                   
        if not cols:
            # I need to generate a series of colours.
            cmap = cm.get_cmap(cm.Paired, len(kargs["cond_names"]))
            cols = []
            step = 256 // len(kargs["cond_names"]) 
            for t in xrange(1, 256, step):
                cols.append(cmap(t))
            #print cols
        
        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.3, 0.1, 0.68, 0.8]) # plenty of space for labels
        
        # Convert to Numpy arrays for type laziness
        da = array(da).T
        err = array(err).T
        err_up = array(err_up).T
        err_dn = array(err_dn).T
        wid = (1.0 / len(da))-0.05
        x = arange(len(da[0]))
               
        if "cond_names" in kargs:
            labs = kargs["cond_names"]
        else: # fake one
            labs = ["" for t in da]
    
        general_args = {"ec": "black", "ecolor": "black"} 
        
        for i, r in enumerate(da):
            if err:
                ax.barh(x+(wid*i), r, wid, xerr=err, label=labs[i], fc=cols[i], **general_args)
            elif "err_up" in kargs and "err_dn" in kargs:
                ax.barh(x+(wid*i), r, wid, xerr=(err_dn[i], err_up[i]), label=labs[i], fc=cols[i], **general_args)
            elif "err_up" in kargs:
                ax.barh(x+(wid*i), r, wid, xerr=err_up[i], label=labs[i], fc=cols[i], **general_args)
            else:
                ax.barh(x+(wid*i), r, wid, label=labs[i], fc=cols[i], **general_args)  
            # I'm sure you don't mean just err_dn
        ax.set_ylim([0, x[-1]+(wid*2)])
        
        if "cond_names" in kargs:
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
            
        if "labels" in kargs and kargs["labels"]:
            ax.set_yticklabels(genelist[kargs["labels"]], rotation="horizontal")
            ax.set_yticks(x+0.5)
        
        self.do_common_args(ax, **kargs)
        
        [item.set_markeredgewidth(0.0) for item in ax.xaxis.get_ticklines()]
        [item.set_markeredgewidth(0.0) for item in ax.yaxis.get_ticklines()]
        #[t.set_fontsize(6) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(6) for t in ax.get_yticklabels()]
        
        return(self.savefigure(fig, filename))

    def pie(self, data, labels, filename, aspect="square", title=None, colours=None, 
        draw_percents=False, cmap=None, **kargs):
        """
        Draw a PIE!
        """
        
        fig = self.getfigure(aspect=aspect, **kargs)
        ax = fig.add_subplot(111)
        
        extra_args = {}
        if draw_percents:
            extra_args["autopct"] = "%1.1f%%"
        
        if colours:
            ax.pie(data, labels=labels, shadow=False, colors=colours, **extra_args)
            
        elif cmap:
            ld = len(data)
            cmap = cm.get_cmap(cmap, ld)
            ran = np.linspace(0.0, 1.0, num=ld)#/float(ld)
            colours = cmap(ran)[:,:-1]
            colours = [rgb2hex(col) for col in colours] # rgb2hex from matplotlib not needed, but easier to print
            
            ax.pie(data, labels=labels, shadow=False, colors=colours, **extra_args)           
        else:
            ax.pie(data, labels=labels, shadow=False, **extra_args)
        
        if title:
            ax.set_title(title)
        
        real_filename = self.savefigure(fig, filename)
        return(real_filename)