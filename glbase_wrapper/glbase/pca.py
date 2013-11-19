"""

PCA analysis for glbase expression objects.

"""

from operator import itemgetter

import numpy
import numpy.linalg as LA
import matplotlib.pyplot as plot
import matplotlib.patches
# not always available?
from mpl_toolkits.mplot3d import Axes3D, art3d

import config
from draw import draw
from genelist import genelist

class pca:
    def __init__(self, parent=None, rowwise=False, label_key=None):
        """
        **Purpose**
            A custom class for PCA analysis of (at the moment) expression object data.
            
            Not recommended to use directly, but if you must then:
            
            expn = expression(...)
            expn_pca = pca(expn)
            
            expn_pca.loading(...)
            ...
            
        **Arguments**
            parent (Required)
                The parent expression/genelist object.
                
            rowwise (Optional, default=False)
                Perform PCA on the columns
                
            label_key (Optional, Required if rowwise=True)
                The key to use in the genelist to label 
        """
        self.parent = parent
        
        matrix = numpy.array(parent.serialisedArrayDataList)
               
        self.__draw = draw()
        self.cols = "black"
        self.rowwise = rowwise
        
        if rowwise:
            self.__u, self.__d, self.__v = LA.svd(matrix, full_matrices=False)
            self.labels = parent[label_key]
        else:
            self.__u, self.__d, self.__v = LA.svd(matrix.T, full_matrices=False)
            self.labels = parent.getConditionNames()
            
        self.__u = numpy.array(self.__u)
        self.valid = True # Just check it's all calc'ed.

    def __repr__(self):
        return("<glbase.pca>")
        
    def __str__(self):
        ret = ["PCA object",
            "    Expression: %s" % self.parent.name,
            "    Dimensions: %s" % len(self.__d),
            "    u: %s" % str(self.__u.shape),
            "    v: %s" % str(self.__v.shape)
            ]
        return("\n".join(ret))

    def __len__(self):
        """
        (Override)
        return the number of dimensions.
        """
        return(len(self.__d))
        
    def get_uvd(self):
        """
        **Purpose**
            Get the u, v, d matrices.
            
        **Arguments**
            None
            
        **Returns**
            A dict, containing:
            {"u": u,
            "v": v,
            "d": d}
        """
        ret = {"u": self.__u,
            "v": self.__v,
            "d": self.__d}
        return(ret)
    
    def max(self):
        """
        **Purpose**
            Return the maximum number of PC for this experiment.
        
        **Arguments**
            None
            
        **Returns**
            The number of PCs. Strictly, len(d)
        """
        return(len(self.__d))

    def set_cols(self, sample_colours):
        """
        **Purpose**
            Set the colours for the individual samples.
            
        **Arguments**
            sample_colours (Required)
                A list of sample colours (or None) for scatter plots and the like.
                Must be the same length as the number of conditions.
                
        **Returns**
            None
        """
        self.cols = sample_colours
    
    def pcloading(self, filename=None, **kargs):
        """
        **Purpose**
            plot a graph of PC loading.
            
        **Arguments**
            filename (Required)
                The filename to save the image to.
                
            <common figure arguments are supported>
        
        **Returns**
            None.
        """
        assert filename, "loading(): Must provide a filename"
        
        if "aspect" not in kargs:
            kargs["aspect"] = "wide"
        
        fig = self.__draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        x = numpy.arange(len(self.__d))
        ax.bar(x-0.4, self.__d, ec="black", color="grey")
        ax.set_xlabel("Principal components")
        ax.set_ylabel("Loading")
        ax.set_xticklabels(x+1)
        ax.set_xticks(x)
        ax.set_xlim([-0.5, len(self.__d)-0.5])
        self.__draw.do_common_args(ax, **kargs)
        real_filename = self.__draw.savefigure(fig, filename)
        config.log.info("loading(): Saved PC loading '%s'" % real_filename)
        
    def scatter(self, x, y, filename=None, spot_cols=None, label=False, alpha=0.8, 
        spot_size=40, spot_label_textsize=7, cut=None, squish_scales=False, **kargs): 
        """
        **Purpose**
            plot a scatter plot of cond1 against cond2.
        
        **Arguments**
            x, y (Required)
                PC dimension to plot as scatter
                Note that PC begin at 1 (and not at zero, as might be expected)
            
            filename (Required)
        
            spot_cols (Optional, default="black" or self.set_cols())
                list of colours for the samples, should be the same length as 
                the number of conditions. 
            
            label (Optional, default=False)
                label each spot with the name of the condition
                
            alpha (Optional, default=0.8)
                alpha value to use to blend the individual points
                
            spot_size (Optional, default=40)
                Size of the spots on the scatter
                
            spot_label_textsize (Optional, default=7)
                Size of the spot label text, only valid if label=True
        
            cut (Optional, default=None)
                Send a rectangle of the form [topleftx, toplefty, bottomrightx, bottomrighty], cut out all of the items within that
                area and return their label and PC score    
                
            squish_scales (Optional, default=False)
                set the limits very aggressively to [min(x), max(x)]
        
        **Returns**
            None
            You can get PC data from pca.get_uvd()
        """
        assert filename, "scatter(): Must provide a filename"     

        ret_data = None      
        xdata = self.__v[x-1]
        ydata = self.__v[y-1]
        
        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        
        fig = self.__draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        cols = self.cols
        if spot_cols:
            cols = spot_cols            
        
        ax.scatter(xdata, ydata, s=spot_size, alpha=alpha, edgecolors="none", c=cols)
        if label:
            for i, lab in enumerate(self.labels):
                ax.text(xdata[i], ydata[i], lab, size=spot_label_textsize, ha="center", va="top")
        
        # Tighten the axis
        if squish_scales:
            if not "xlims" in kargs:
                ax.set_xlim([min(xdata), max(xdata)])
        
            if not "ylims" in kargs:
                ax.set_ylim([min(ydata), max(ydata)])
        
        ax.set_xlabel("PC%s" % (x,)) # can be overridden via do_common_args()
        ax.set_ylabel("PC%s" % (y,))
        
        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])
        
        if cut:
            rect = matplotlib.patches.Rectangle(cut[0:2], cut[2]-cut[0], cut[3]-cut[1], ec="none", alpha=0.2, fc="orange")
            ax.add_patch(rect)

            tdata = []
            for i in xrange(0, len(xdata)):
                if xdata[i] > cut[0] and xdata[i] < cut[2]:
                    if ydata[i] < cut[1] and ydata[i] > cut[3]:
                        tdata.append({"name": self.labels[i], "pcx": xdata[i], "pcy": ydata[i]})
            if tdata:
                ret_data = genelist()
                ret_data.load_list(tdata)
            
        self.__draw.do_common_args(ax, **kargs)
        
        real_filename = self.__draw.savefigure(fig, filename)
        config.log.info("scatter(): Saved 'PC%s' vs 'PC%s' scatter to '%s'" % (x, y, real_filename)) 
        return(ret_data)

    def scatter3d(self, x, y, z, filename=None, spot_cols=None, label=False, stem=False, 
        label_font_size=6, rotation=134, interactive=False, **kargs): 
        """
        **Purpose**
            plot a scatter plot of PC1 against PC2 against PC3
            This is the 3D variant
                    
        **Arguments**
            x, y, z (Required)
                PC dimension to plot as scatter
                Note that PC begin at 1 (and not at zero, as might be expected)
        
            spot_cols (Optional, default="black" or self.set_cols())
                list of colours for the samples, should be the same length as 
                the number of conditions. 
            
            label (Optional, default=False)
                label each spot with the name of the condition
            
            label_font_size (Optional, default=6)
                label font size.
            
            stem (Optional, default=False)
                Draw stems from the point down to the base of the graph. (i.e. z=0)
                
            rotation (Optional, default=134)
                x,y, plane rotation
                
            interactive (Optional, default=False)
                if True then spawn the matplotlib show() view. Note that this will pause
                execution of your script.
                Note that by default glbase uses a non-GUI matplotlib setting.
                
                You will need to fiddle around with matplotlib.use() before importing glbase         
        
        **Returns**
            None
            You can get PC data from pca.get_uvd()
        """
        assert filename, "scatter(): Must provide a filename"     
      
        xdata = self.__v[x-1]
        ydata = self.__v[y-1]
        zdata = self.__v[z-1]
            
        fig = self.__draw.getfigure(**kargs)
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=rotation)
        
        cols = self.cols
        if spot_cols:
            cols = spot_cols            
        
        ax.scatter(xdata, ydata, zdata, edgecolors="none", c=cols, s=40)
        if label:
            for i, lab in enumerate(self.labels):
                ax.text(xdata[i], ydata[i], zdata[i], lab, size=label_font_size, ha="center", va="bottom")
        
        if stem: # stem must go after scatter for sorting
            z_min = min(zdata)
            for x_, y_, z_ in zip(xdata, ydata, zdata):        
                line = art3d.Line3D(*zip((x_, y_, z_min), (x_, y_, z_)), marker=None, c="grey", alpha=0.5)
                ax.add_line(line)
        
        ax.set_xlabel("PC%s" % (x,)) # can be overridden via do_common_args()
        ax.set_ylabel("PC%s" % (y,))
        ax.set_zlabel("PC%s" % (z,))
        
        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])
        
        self.__draw.do_common_args(ax, **kargs)
        
        if interactive:
            fig.show() # hope you are not on a cluster!
                
        real_filename = self.__draw.savefigure(fig, filename)
        
        config.log.info("scatter3d(): Saved 'PC%s' vs 'PC%s' vs 'PC%s' scatter to '%s'" % (x, y, z, real_filename))     

    def gene_loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, **kargs):
        """
        Deprecated.
        
        NOTE: This is an alias of loading(), please use loading() and NOT gene_loading()
        in future this method will be deleted.
        """
        return(self.loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, **kargs))

    def loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=None, **kargs):
        """
        **Purpose**
            Get the loading for the items for a particular PC
            
            (technically, get u for a particular PC)
            
        **Arguments**
            filename(Optional)
                filename to save the loading barchart to. This is optional, so you can use this function
                just top return top and bot
                
            PC (Required)
                The Principal Component to use
                
            label_key (Required)
                The key in the expression object to use to label the bar chart
            
            top & bot (Optional, default=50)
                the number of top and bottom genes to plot on the bar graph and
                also to return as a genelist. Set both to None to get all componenets
                
            all (Optional, default=False)
                if all is True, return all of the items.
        
        **Returns**
            topbot of the loading in a new genelist with an extra key "loadingPC<PC#>"
            
            The object will be based on the original expression object, and will be sorted
            based on the PC loading. This means you can do some neat stuff like:
            
            new_expn = pca.gene_loading(..., top=50, bot=50)
            new_expn.heatmap()
            new_expn.boxplot()
            
        """
        if not self.rowwise:
            return(self.row_loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, all=all, **kargs))
        else:
            return(self.condition_loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, all=all, **kargs))
        
    def row_loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=False, **kargs):
        """
        Deprecated, please use loading()
        """
        assert not self.rowwise, "row_loading(): You probably don't mean to use this when rowwise=True, try condition_loading()"
        assert PC >= 1, "row_loading(): PC of <1 specified"
        
        if "aspect" not in kargs:
            kargs["aspect"] = "long"
        
        data = self.__u[:,PC]
        labs = self.parent[label_key]
        packed_data = [{label_key: i[0], "l": i[1]} for i in zip(labs, data)]
        
        sorted_data = sorted(packed_data, key=itemgetter("l"))
        data = [i["l"] for i in sorted_data]
        labs = [i[label_key] for i in sorted_data]
        
        if all:
            data = data
            labs = labs
        else:
            if bot > 0 and top > 0: # data[-0:] returns the entire list and data[0:0] returns [] !
                data = data[0:top] + data[-bot:]
                labs = labs[0:top] + labs[-bot:]
            elif top > 0:
                data = data[0:top]
                labs = labs[0:top]             
            elif bot > 0:
                data = data[-bot:]
                labs = labs[-bot:]        
        
        if filename:
            fig = self.__draw.getfigure(**kargs)
            ax = fig.add_subplot(111)
            ax.set_position([0.3,0.03,0.6,0.96])
        
            x = numpy.arange(len(data))
            ax.barh(x-0.4, data, ec="black", color="grey")
            ax.set_ylabel("Rows")
            ax.set_xlabel("Loading")
            ax.set_yticklabels(labs)
            ax.set_yticks(x)
            ax.set_ylim([-0.5, len(data)-0.5])
            [t.set_fontsize(6) for t in ax.get_yticklabels()]
        
            self.__draw.do_common_args(ax, **kargs)
            real_filename = self.__draw.savefigure(fig, filename)
        
            config.log.info("row_loading(): Saved PC gene_loading '%s'" % real_filename)
        
        # work out the list to return
        newgl = genelist()
        newgl.load_list([{label_key: i[0], "pc_loading": i[1]} for i in zip(labs, data)]) # relist it so that top bot are used
        newexpn = newgl.map(genelist=self.parent, key=label_key, greedy=False)
        newexpn.sort("pc_loading")
        return(newexpn)

    def condition_loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=False, **kargs):
        """
        deprecated, please use loading()
        """
        assert self.rowwise, "condition_loading(): You probably don't mean to use this when rowwise=False, try condition_loading()"
        assert PC >= 1, "condition_loading(): PC of <1 specified"
        
        if "aspect" not in kargs:
            kargs["aspect"] = "long"
        
        data = self.__u[:,PC]
        labs = self.parent._conditions
        packed_data = [{label_key: i[0], "l": i[1]} for i in zip(labs, data)]
        
        sorted_data = sorted(packed_data, key=itemgetter("l"))
        data = [i["l"] for i in sorted_data]
        labs = [i[label_key] for i in sorted_data]
        
        if all:
            data = data
            labs = labs
        else:
            if bot > 0 and top > 0: # data[-0:] returns the entire list and data[0:0] returns [] !
                data = data[0:top] + data[-bot:]
                labs = labs[0:top] + labs[-bot:]
            elif top > 0:
                data = data[0:top]
                labs = labs[0:top]             
            elif bot > 0:
                data = data[-bot:]
                labs = labs[-bot:]        
        
        if filename:
            fig = self.__draw.getfigure(**kargs)
            ax = fig.add_subplot(111)
            ax.set_position([0.3,0.03,0.6,0.96])
        
            x = numpy.arange(len(data))
            ax.barh(x-0.4, data, ec="black", color="grey")
            ax.set_ylabel("Columns")
            ax.set_xlabel("Loading")
            ax.set_yticklabels(labs)
            ax.set_yticks(x)
            ax.set_ylim([-0.5, len(data)-0.5])
            [t.set_fontsize(6) for t in ax.get_yticklabels()]
        
            self.__draw.do_common_args(ax, **kargs)
            real_filename = self.__draw.savefigure(fig, filename)
        
            config.log.info("condition_loading(): Saved PC gene_loading '%s'" % real_filename)
        
        # work out the list to return
        newgl = genelist()
        newgl.load_list([{"name": i[0], "pc_loading": i[1]} for i in zip(labs, data)]) # relist it so that top bot are used
        newgl.sort("pc_loading")
        return(newgl)
        