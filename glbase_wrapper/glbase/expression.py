"""
**Purpose**
    An all-purpose container for transcript expression data.

**to do**

* conditions are not always stored as floats... [fixed?]
* need to abstract away to draw.
* often the reader will choke on an array... [partially fixed?]
* move the "draw" stuff into a seperate heritable class. [partially implemented]
* to draw, do things like: array.draw.heatmap()
* normaliseToCondition has a bug which is currently bodged over - if the background = 0.0 then it will DivByZero Error.
* a 'expression' class that contains useful methods to deal with microarray data, rather than the
    very thin implementation here of a list buried in conditions
    (surely a numpy array would be more approporiate?)
* scipy.cluster relies of recursion to perform the clustering.
    On large sized arrays it chokes horribly.
"""

import sys, os, csv, string, math, copy

import numpy
from numpy import array, arange, meshgrid, zeros, linspace, mean, object_, std
from array import array as qarray # this should be deprecated later.

import config
from flags import *
from base_expression import base_expression
from draw import draw
from progress import progressbar
from errors import AssertionError, ArgumentError
from genelist import genelist
from location import location

import matplotlib.pyplot as plot
from pylab import bivariate_normal, griddata # comes from where?
import matplotlib.cm as cm

class expression(base_expression):
    """
    **Purpose**

        A container class for transcript expression data. For example, microarrays, RNA-seq
        or any data with measurable transcript abundance.
        
        You should be aware that expression data is very basic in glbase. And 
        expression-data analysis in glbase is very underpowered and requires normalised
        and easily manipulatable data. Examples for normalisation are
        genespring output, or output from R and PMA. For RNA-seq, cufflinks output, etc.
        
        Useful functions for manipulating expression data in
        relation to other genelists.
        
        Inherits all genelist methods and implements several expression
        specific methods.

    **Arguments**
        filename (Required, one of loadable_list or filename)
            the filename of the expression-data to load.

        loadable_list (Required, one of loadable_list or filename)
            a genelist-like object I can use to construct the list from. If you use this then
            expn should be a list of keys to extract the expression data from.
            
        format (Required)
            a format specifier.
            Must include some form of {"conditions": {"code": "column[4:]"}}
            to specifiy the location of the numeric array data.
            glbase will try to guess but will likely fail.

            There are a few special import methods available.
            To import illumina microarray data output from bead studio
            use the special format specifier: "illumina"

        expn (Required)     
        
        	If filename is being used:
            Some sort of descriptor telling me where the expression data actually is.
            For example:
                "column[3:]" (Take each column from column3 onwards, until the end column)
            Or:
                "column[4::2]" (Take every second item, from column 4 until the end)
            Or:
                "column[5:8]" (Take column 5 through 7 - 8 is not a typo. The lists are 
                    zero-ordered and closed)
                    
            In fact you can even give compound statements:
                "[column[7], column[17]]" (Take columns 7 and 17)
                
            The only rules are it must be a valid piece of python code.
            
        	If a loadable_list:
            This should be the name of the keys used to extract the expresion data
    """
    def __init__(self, loadable_list=None, filename=None, format=None, expn=None, **kargs):
        if loadable_list:
            base_expression.__init__(self, loadable_list=loadable_list, expn=expn)
        else:
            base_expression.__init__(self, filename=filename, expn=expn, format=format, **kargs)

    def __repr__(self):
        return("glbase.expression")

    def findGene(self, **kargs):
        """
        # valid keys are in do
        accepts things like:
        findGene(refseq="123344", symbol="Pou5f1")
        (Obselete?)
        """
        ret = []
        do = frozenset(["refseq", "entrez", "symbol", "coords"])

        for key, value in kargs.iteritems():
            if key == "loc":
                ret.append(self._findByCoords(key, value))
            else:
                try:
                    ret.append(self._findDataByKeyGreedy(key, value))
                except KeyError:
                    print "Warning: key is not present in array"
        return(ret)

    def getConditionNames(self):
        """
        returns a list of the condition headers
        """
        return(self.conditions)

    # ----------- overrides/extensions ---------------------------------

    def getArrayName(self):
        return(self.array["name"])

    def getGenomeName(self):
        if not self.genome:
            return("Genome not bound")
        return(self.genome.getName())

    def sliceConditions(self, conditions=None, **kargs):
        """
        **Purpose**

            return a copy of the expression-data, but only containing
            the condition names specified in conditions

        **Arguments**

            conditions (Required)

                A list, or other iterable of condition names to extract
                from the expression-data. Every condition name must be present
                on the expression-data.

        **Result**

            A new expression-data object with the same settings as the original,
            but containing only the expression-data conditions specified in
            the 'conditions' argument.
        """
        valig_args = ["conditions"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.sliceConditions, k)

        assert conditions, "You must specify a list of conditions to keep"
        assert isinstance(conditions, list), "You must specify a list of conditions to keep"
        for item in conditions:
            assert item in self.conditions, "'%s' condition not found in this expression-data" % item

        newl = copy.deepcopy(self)

        copymask = [] # work out a mask to extract the correct array columns.
        for name in conditions:
            for i, c in enumerate(self.conditions):
                if c == name:
                    copymask.append(i)

        for index, item in enumerate(self.linearData):
            new_array_data = []
            for c in copymask: # use the mask to collect the array entries.
                new_array_data.append(item["conditions"][c])
            newl.linearData[index]["conditions"] = new_array_data

        newl._history.append("sliced conditions, kept %s" % "".join(conditions))
        newl.conditions = conditions
        newl._optimiseData()
        return(newl)

    def getDataForCondition(self, condition_name):
        """
        **Purposse**
            get all of the expression-data data for a particular condition
            name, returns a list of all the values.
            The list remains in the same order as the overall list,
            so to get e.g. the refseq column, do: l = microarray["refseq"]
        """
        names = self.getConditionNames()
        assert condition_name in names, "No condition named '%s' on this array" % condition_name

        l = []
        con = self["conditions"] # get all of the array data
        names = self.getConditionNames()
        index = names.index(condition_name)

        return([i[index] for i in con])

    def drawHeatmap(self, filename=None, row_label_key="name", **kargs):
        """
        **Purpose**

        draw a simple heatmap of the current expression-data data.

        **Arguments**

        filename (Required)
            the filename of the image to save. depending upon the current
            drawing settings it will save either a png (default) svg or eps.

        bracket (Optional, default = no bracketing performed)
            bracket the data within a certain range of values.
            For example to bracket to 0 .. 1 you would use the syntax::

                result = array.drawHeatmap(filename="ma.png", bracket=[0,1])

            Or for something like log2 normalised array data:

                result = array.drawHeatmap(filename="ma.png", bracket=[-2,2])

            "Bracket' chops off the edges of the data, using this logic:

                if value > high_bracket then value := high_bracket
                if value < low_bracket then value := low_bracket

            See normal for a method that modifies the data.

        row_label_key (Optional, default="name")
            A key in your genelist to use to label the rows. Examples would be gene names accesion
            numbers or something else.

        normal (Optional, default = no normalising)
            Unimplemented

        row_cluster (Optional, default = True)
            cluster the rows? True or False

        col_cluster (Optional, default = True)
            cluster the column conditions, True or False

        log (Optional, defualt=False, True|False of 2..n for log2, log10)
                log the y axis (defaults to e)
                send an integer for the base, e.g. for log10

                log=10

                for log2

                log=2

                for mathematical constant e

                log=True
                log="e"

        **Result**

        saves an image to the 'filename' location and
        returns the 'actual filename' that really gets saved (glbase
        will modify e.g. a '.png' ending to '.svg' or '.eps' etc. depending
        upon the current setings).
        """
        # checks for option here please.
        assert filename, "you must specify a filename"

        data = self.serialisedArrayDataList
        if "log" in kargs:
            data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])

        # convert it into the serialisedArrayDataDict that _heatmap() expects.
        newdata = {}
        for index, name in enumerate(self.conditions):
            newdata[name] = data[index] # get the particular column

        actual_filename = self.draw.heatmap(data=newdata,
            row_names=self[row_label_key], col_names=self.getConditionNames(),
            filename=filename, **kargs)

        config.log.info("Saved the heatmap image: %s" % actual_filename)
        return(actual_filename)

    def normaliseToCondition(self, condition_name, bUseFoldChange=True, keep_normed=False, **kargs):
        """
        **Purpose**
            normalise all other conditions to condition_name and delete
            condition name from the expression-data list

        **Arguments**
            keep_normed (boolean, Optional)
                keep the original data, but set it to 0

            ??? not true right:
            This puts the data in the range 0 .. 1 where 0.5 = the same as the control.

        **Returns**
            returns the newly normalised list
        """
        names = self.getConditionNames()

        assert condition_name in names, "condition name: %s is not in this expression-data" % condition_name

        name_index = names.index(condition_name)
        #print name_index

        newl = self.__copy__()
        newl.linearData = []
        newl.conditions = []

        p = progressbar(len(self.linearData))
        for index, item in enumerate(self.linearData):
            #print item
            old_array_data = item["conditions"]
            #print old_array_data
            new_array_data = []
            toNormal = old_array_data[name_index]+0.0000001 # stop divByZero errors.
            for i, datum in enumerate(old_array_data):
                if names[i] != condition_name:
                    if bUseFoldChange:
                        new_array_data.append(math.pow(2, -(float(toNormal) / (float(datum)+0.0000001))))
                    else:
                        new_array_data.append(float(toNormal) / (float(datum)+0.0000001))
                elif keep_normed:
                    new_array_data.append(1.0) # er... this should be 0.5 ?

            data_copy = copy.deepcopy(item)
            newl.linearData.append(data_copy)
            data_copy["conditions"] = new_array_data # load the new_array_data over the old one.

            p.update(index)

        # rebuild the condition names (optimiseData can't handle this)
        if keep_normed:
            newl.conditions = names
        else:
            # delete the old label
            newl.conditions = []
            for name in names:
                if name != condition_name:
                    newl.conditions.append(name)
                elif keep_normed:
                    newl.conditions.append(name)

        newl._optimiseData()
        config.log.info("Normalised to condition: %s" % condition_name)
        self._history.append("Normalised To Condition: %s" % condition_name)
        return(newl)

    def drawDotplot(self, x_condition_name, y_condition_name, filename=None, genelist=None, key=None, **kargs):
        """
        **Purpose**
            draw an X/Y dot plot or scatter plot, get R^2 pearson correlation etc.

        **Arguments**
            x_condition_name = the name of the er... X condition
            y_condition_name = the name of the er... Y condition
            
            genelist (Optional)
                If you send a peaklist and a key then these items will be emphasised on the dotplot.
                
            key (Optional, Required if 'genelist' used)
                The key to match between the expression data and the genelist.
                
            available key-word arguments:
            xlabel, ylabel, title, log (set this to the base to log the data by),
            xlims, ylims

        **Returns**
            the actual filename saved as and a new image in filename.
        """
        assert filename, "no filename specified"
        assert x_condition_name in self.serialisedArrayDataDict, "%s x-axis condition not found" % x_condition_name
        assert y_condition_name in self.serialisedArrayDataDict, "%s y-axis condition not found" % y_condition_name

        x_data = self.getDataForCondition(x_condition_name)
        y_data = self.getDataForCondition(y_condition_name)
               
        if genelist and key:
            matches = genelist.map(genelist=self, key=key) # make sure resulting object is array
            tx = matches.getDataForCondition(x_condition_name)
            ty = matches.getDataForCondition(y_condition_name)
            
            real_filename = self.draw.nice_scatter(x_data, y_data, filename, spots=(tx, ty), **kargs)
        else:
            real_filename = self.draw.nice_scatter(x_data, y_data, filename, **kargs)
            
        """
        if len(x_data) < 100:
            # for the matplot.lib < 100: I want to label everything.
            names = self["name"]
            for i, n in enumerate(names):
                plot.annotate(n, (x_data[i], y_data[i]), size=6)
        """

        config.log.info("Saved the scatter plot as '%s'" % real_filename)
        return(True)

    def drawBoxplot(self, filename=None, **kargs):
        """
        **Purpose**

        Draw a boxplot of all conditions.

        **Arguments**

            filename (Required)
                filename to save as. The file extension may be modified
                depending the setting of the current config.DEFAULT_DRAWER

            log (True|False of 2..n for log2, log10)
                log the y axis (defaults to e)
                send an integer for the base, e.g. for log10

                log=10

                for log2

                log=2

                for mathematical constant e

                log=True
                log="e"

        **Results**

        saves an image with the correct filetype extension for the current
        config.DEFAULT_DRAWER.
        returns the actual filename used to save the file.
        """
        assert filename

        ylimits = False
        do_log = False

        data = self.serialisedArrayDataList
        for key in kargs:
            if key == "log":
                data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])
            if key == "ylimits": # requires a tuple.
                ylimits = kargs[key]
            if key == "filename":
                filename = kargs[key]

        # do plot
        actual_filename = self.draw.boxplot(data=data, filename=filename, labels=self.getConditionNames())

        config.log.info("Saved the boxplot image: %s" % actual_filename)
        return(actual_filename)

    def __log_transform_data(self, serialisedArrayDataList=None, log=math.e):
        """
        (Internal)

        transforms the data based on base
        helper for drawBoxPlot() and draw drawCurves()
        
        Zeros are trimmed from the data. This is important because now the output
        may not be synchronised to the input. Care should be taken with this transform.
        """
        assert serialisedArrayDataList, "[Internal] __log_transform_data() - no data provided"

        do_log = False

        if log == math.e:
            do_log = math.e
        elif isinstance(log, bool):
            do_log = math.e
        elif isinstance(log, int):
            do_log = log
        else:
            do_log = False

        if do_log:
            data = []
            for set in serialisedArrayDataList:
                row = []
                for index, item in enumerate(set):
                    if item != 0.0:
                        row.append(math.log(item, do_log))
                    else:
                        row.append(math.log(0.000001, do_log)) # Append a very small value
                data.append(row)
            return(data)
        else:
            return(serialisedArrayDataList)

    def drawCurves(self, filename=None, **kargs):
        """
        **Purpose**

        draw a bell-curve diagram of the expression-data expression.

        **Arguments**

            filename (Required)
                filename of the resulting

            window
                size of window for moving average

            modifier
                undocumented fudge for float based arrays.

            xlimits
                a tuple of the form: (minimum_x, maximum_y)

            log (True|False of 2..n for log2, log10)
                log the y axis (defaults to e)
                send an integer for the base, e.g. for log10

                log=10

                for log2

                log=2

                for mathematical constant e

                log=True
                log="e"
                
                Data points that are 0.0 will be trimmed from the data.
                This means the number of samples plotted may not
                be the same as the number of points in the microarray
                data.

            cumulative (True|False, default False)

                draw cumulative curves.
                
            verbose (True|False), default=False
                print out the means and standard deviations.

        **Result**

        saves an image to 'filename'
        returns the actual filename (the filename may be modified
        depending upon the current display driver)

        """
        assert filename, "no filename given"

        window_size = 200
        simple_args = ["filename", "window", "modifier"]
        xlimits = None
        do_log = False
        extra_args = {}
        data = self.serialisedArrayDataList

        for k in kargs:
            if k == "log":
                data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])
            if k == "filename":
                filename = kargs[k]
            if k == "window":
                window_size = kargs[k]
            if k == "modifier": # undocumented fudge for flat arrays.
                binner_modifier = kargs[k]
            if k == "xlimits": # requires a tuple.
                xlimits = kargs[k]
            if k == "cumulative":
                extra_args["cumulative"] = kargs["cumulative"]

        # normalise data, for each condition.
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        if "verbose" in kargs and kargs["verbose"]:
            print "name\tmean\tstd" 

        for i, c in enumerate(self.getConditionNames()):
            ax.hist(data[i], bins=window_size, histtype="step", label=c, **extra_args)
            m = mean(data[i])
            d = std(data[i])
            ax.axvline(x=m, color="red")
            ax.axvline(x=m-d, color='grey', ls=":")
            ax.axvline(x=m+d, color='grey', ls=":")
            if "verbose" in kargs and kargs["verbose"]:
                print "%s\t%.2f\t%.2f" % (c, m, d)

        if xlimits: 
            ax.xlim(xlimits)
        ax.legend()
        
        real_filename = self.draw.savefigure(fig, filename)
        config.log.info("Saved the Distribustion Curves: %s" % filename)
        return(True)

    def getDataByCriteria(self, **kargs):
        """
        **Purpose**
            used to test for fold_up, sig_up, etc...
    
            function can be any helper_function, which uses data[] as it's set for each item.
            a set of already defined functions exist in flags.py
            function must accept these arguments: (data[], conditionNames)
            data is a dict of the form {"con_name1": value, "con_name2": value ...}
        """
        function = None
        normal = None
        for k in kargs:
            if k == "function":
                function = kargs[k]
            if k == "normal":
                normal = kargs[k]

        if not function:
            config.log.error("Criteria function unavailable.")
            return(False)

        newl = self.__copy__()
        newl.linearData = []

        conNames = self.getConditionNames()
        for item in self:
            data = item["conditions"]
            # package as a dict:
            dd = {}
            for index, name in enumerate(conNames):
                dd[name] = data[index]
            if function(dd, conNames, normal, **kargs): # pass on other kargs
                newl.linearData.append(deepcopy(item))

        newl._optimiseData()
        return(newl)

    def _insertCondition(self, condition_name, condition_data, range_bind=None, **kargs):
        """
        (Internal)
        candidate for reveal?
        """
        self.conditions.append(condition_name)
        max_data = max(condition_data)
        min_data = min(condition_data)
        if len(condition_data) != len(self):
            config.log.error("Insertion of array data, wrongly sized")
            return(False)
        for index, item in enumerate(self):
            if range_bind:
                toAdd = (float((condition_data[index] - min_data)) / max_data)
                toAdd = (toAdd * (range_bind[0] + range_bind[1])) - range_bind[0]
            else:
                toAdd = condition_data[index]
            item["conditions"].append(toAdd)
        self._optimiseData()
        return(True)

    def remove_by_expression(self, minimum_value, number_to_pass=1, **kargs):
        """
        **Purpose**
            If the probe expression is below 'minimum_value' in at least
            'number_to_pass' conditions

        **Arguments**

        **Returns**
            A new expression object with probes that fail to pass removed.
        """
        # kargs tidier to go here.

        newl = self.__copy__()
        newl.linearData = []

        for probe in self.linearData:
            s = sum([i > minimum_value for i in probe["conditions"]])
            if s > number_to_pass:
                newl.append(probe)
        newl._optimiseData()
        config.log.info("'%s' probes failed the criteria" % (len(self) - len(newl)))
        newl._history.append("'%s' probes failed the criteria" % (len(self) - len(newl)))
        return(newl)

    def cumulative_distributions(self, genelists=None, filename=None, key=None, **kargs):
        """
        **Purpose**

            draw a set of cumulative distributions, based on a selection
            of genelist-like objects that can be mapped to
            this microarray using 'key'

        **Arguments**

            genelists (Required)
                a list or other iterable of genelists

            filename (Required)
                the filename to save the image to.

            key (Required)
                the key to use to match the microarray to the genelist.

        **Returns**

            An image, saved to 'filename'
        """
        valig_args = ["genelists", "filename", "key"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.cumulative_distributions, k)

        assert filename, "you must specify a valid filename"
        assert key, "you must specify a mapping key"
        assert genelists[0], "you must specify a valid list of genelists"
        assert key in genelists[0], "key '%s' not found in the genelists" % key
        assert key in self, "key '%s' not found in the expression-data" % key

        mapped_scores = []

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        axis = fig.add_subplot(111)

        for g in genelists:
            mapped = self.map(genelist=g, key=key)

            nmap = []
            # this will sum all items in array.
            for a in mapped:
                print a["conditions"]
                s = sum(a["conditions"])
                nmap.append(s)

            # cumulate the nmap
            for i, v in enumerate(nmap):
                try:
                    nmap[i] = nmap[i] + nmap[i+1]
                except:
                    break

            print nmap

        # matplotlib junk is inappropriately here: to go later.

        axis.plot(nmap, label=g.name)

        axis.set_title("")
        #axis.show_legend()
        fig.savefig(filename)

    def drawBarChart(self, gene_symbols=None, filename=None, key=None, labels=None,
        errs_are_absolute=False, error_keys=None, fake_entries=True, **kargs):
        """
        **Purpose**
            Draw barcharts of the individual genes or a (more likely) list of genes.
            I will use these symbols and collect every matching entry in the expression-data 
            and will draw a horizontal bar-chart
            
        **Arguments**
            gene_symbols (Required)
                a list or single item that I can find in the expression data.
                
            filename (Required)
                the filename to save the image to.
                
            key (Optional, default="name"|"symbol")
                The key to match the gene_symbols in. 
                Ususally you will mean either name, gene symbol, or enst or equivalent
                annotation
            
            labels (Optional, default=<the same as 'key' argument, above>)
                specify an optional key to use for the actual labels on the chart
                If left unspecified then the 'key' argument is used for the labels.
                This allows you to search using 'key' but actually label the items with 'label'
            
            error_keys (Optional, default=None)
                If the array data has error values then you can send the key names.
                The input comes in two forms:
                
                error_keys="stderr" # this will use one key that should match the 
                                    # expression data
                
                error_keys=["conf_lo", "conf_hi"] # If you have paired keys for 
                    # error values. In this instace conf_lo holds the lower bound confidence
                    # interval and conf_hi the Higher bound.
            
            errs_are_absolute (Optional, default=False)
                If your error bars are absolute then set this to True.
                What it means is if your expression is (say) 10, and your error is 2.
                then the error bars will be drawn from 8 -> 12 (ie. the errors are relative
                to the expression value. If however you are using 
                
            fake_entries (Optional, default=True)
                If the value you search for (usually a gene name) does not exist in the data
                then bar_chart will fake an empty entry, and the gene will still show up on the
                graph even though there is no data in your list. It will be given expression values of 
                0 (actually, a very small float to avoid log errors).
                Set this to False if you don't want empty values to be placed on the bar_chart              
        
            Other arguments understood by the figure:
                title, xlabel, ylabel
                
                and 'cols' - where you can send a list of colours for the bars. The list
                must be as long as the data, or it will throw an error.
        
        **Example**
            Okay, this is pretty complicated (It's because the expression system has only a rudimentary 
            understanding of error values, I may change this in the future)...
            
            x = expression(...)
            
            # I don't have any error values:
            x.bar_chart(filename="meh1.png", gene_symbols=["Nanog", "Sox2", "Sox17", "Stat3"],
                key="name")
            
            # My errors are confidence intervals relative to the expression value
            x.bar_chart(filename="meh2.png", gene_symbols=["Nanog", "Sox2", "Sox17", "Stat3"],
                key="name", error_keys=["conf_lo", "conf_hi"], errors_are_absolute=True)
                
            # My errors are standard errors (or equivalent) and are symmetric around the mean
            x.bar_chart(filename="meh3.png", gene_symbols=["Nanog", "Sox2", "Sox17", "Stat3"],
                key="name", error_keys="stderr")
        
        **Returns**
            A saved image in filename and None
        """
        assert filename, "no filename specified"
        assert gene_symbols, "no gene_symbols specified"
        if "cols" in kargs:
            assert len(kargs["cols"]) == len(self.linearData[0]["conditions"]), "the colour array is not the same length as the array data!"
        
        if not key:
            if "name" in self:
                key = "name"
            else:
                if not "symbol" in self:
                    key = "symbol"
                else:
                    raise AssertionError, "No suitable key found"
        
        if not labels:
            labels = key
        
        if not isinstance(gene_symbols, list):
            gene_symbols = [gene_symbols]
            
        empty_array = [0.0001 for x in self.conditions]
        fake_entry = {"conditions": empty_array, "loc": location(loc="chr1:1000-1000")}
        # work out the error keys:
        if error_keys:
            if isinstance(error_keys, list): # load kargs for draw.bar_chart()
                kargs["err_up"] = error_keys[0]
                kargs["err_dn"] = error_keys[1]
                fake_entry = {"conditions": empty_array, error_keys[0]: empty_array, error_keys[1]: empty_array, "loc": location(loc="chr1:1000-1000")}
            else: # assume string
                kargs["err"] = error_keys
                fake_entry = {"conditions": empty_array, error_keys: empty_array, "loc": location(loc="chr1:1000-1000")}
            
        subset = []
        for g in gene_symbols:
            nl = self._findDataByKeyGreedy(key, g)
            if nl:
                subset = subset + nl
            else:
                new = fake_entry.copy()
                new.update({key: "%s (not detected)" % g, labels: "%s (not detected)" % g}) # I want to signify that these are empty
                subset.append(new)
                config.log.warning("%s not found. Spoofing an empty entry" % g)
        
        if subset:
            nl = genelist()
            nl.load_list(subset)
            
            # This is a kind of peculiar way of doing it...
            self.draw.bar_chart(filename=filename, genelist=nl, data="conditions", 
                errs_are_absolute=errs_are_absolute, 
                labels=labels, cond_names=self.conditions,
                **kargs) 
        
    def hist(self, filename=None, range=None, suppress_zeros=True, log=None, kde=False, 
        covariance=0.2, **kargs):
        """
        **Purpose**
            Draw a normal histogram of the expression values
        
        **Arguments**
            filename (Required)
                the filename to save the resulting image to.
                
            covariance (Optional, default=0.2)
                undocumented wierdness (Sorry)
                
            range (Optional, default=(0, max(expression)))
                This is the maximum value to build the KDE over. (ie. 0 ... mmax).
                You probably really want to chage this value, as small outliers will distort the
                distribution massively.
                Play around with the value until it gets you a nice normal-like distribution.
                
            suppress_zeros (Optional, default=True)
                expression-data are pretty neat and don't include true 0 values for expression. 
                However, RNA-seq must assuredly does include 0 values for some transcripts.
                In some instances you may want to suppress these values.
                
                IF this is set to True expression values of zero are not included in the 
                histogram.
                
            log (Optional, default=False)
                log transform the data.
                At the moment only log2 is supported.
        
            kde (Optional, default=False)
                use kernel density estimation to smooth the data
                
            covariance (Optional, default=0.2)
                Value for KDE smoothing.
            
        **Returns**

        """
        assert filename, "Need filename, leh!"
        
        binned_data = {}
        
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        all_data = {}
            
        for k in self.conditions:
            expn_values = self.getDataForCondition(k)
            if not range: # sample a range if none specified
                range = (0, max(expn_values))
            
            if suppress_zeros:
                expn_values = [v for v in expn_values if int(v*10000) > 0]
            
            expn_values = numpy.array(expn_values)
            
            if log:
                expn_values = numpy.log2(expn_values)
            
            if kde:
                expn_values = utils.kde(expn_values, range=range, covariance=covariance, bins=100)
            
            ax.hist(expn_values, bins=200, range=range, normed=True, histtype='step', label=k)
            all_data[k] = expn_values
            
        ax.legend(ncol=1)
                    
        real_filename = self.draw.savefigure(fig, filename)
        
        config.log.info("Saved '%s'" % real_filename)
        return({"data": all_data, "labels": self.conditions})
        
    def print_means(self, log=None):
        """
        **Purpose**
            Print out the min, max, mean, median and stddev for each condition in the microarray data
        
            TODO: Add test for normality
        
        **Arguments**
            log (Optional)
                if true log2 transform the data.
                
                (Data with an expression of zero will be excluded from the log calculation and
                mean scores)
            
        **Returns**
            a dictionary containing min, max, mean, median, stdev for each condition. In the order
            of "conditions".
        """
        means = []
        medians = []
        stds = []
        mins = []
        maxs = []
        
        for k in self.conditions:
            expn_values = self.getDataForCondition(k)
            if log:
                expn_values = [v for v in expn_values if int(v*10000) > 0]
                expn_values = numpy.log2(expn_values)
            
            aaverage = numpy.average(expn_values)
            amedian = numpy.median(expn_values)
            astd = numpy.std(expn_values)
            amin = self.min()
            amax = self.max()
            
            print "Condition: %s" % k
            print "\tmin", amin
            print "\tmax", amax
            print "\tmean", aaverage
            print "\tmedian", amedian
            print "\tstddev", astd
            
            means.append(aaverage)
            medians.append(amedian)
            stds.append(astd)
            mins.append(amin)
            maxs.append(amax)
            
        return({"mean": means, "median": medians, "stdev": stds, "min": mins, "max": maxs, 
            "conditions": self.conditions})
        
    def min(self):
        """
        **Purpose**
            get the minimum value in the expression data
            
        """
        data = [item["conditions"] for item in self.linearData]       
        data = numpy.array(data)
        return(data.min())
      
    def max(self):
        """
        **Purpose**
            get the maximum value in the expression data
            
        """
        data = [item["conditions"] for item in self.linearData]       
        data = numpy.array(data)
        return(data.max())

    def mean(self):
        """
        **Purpose**
            get the mean value in the expression data
            
        """
        data = [item["conditions"] for item in self.linearData]       
        data = numpy.array(data)
        return(numpy.average(data))
   
    def stdev(self):
        """
        **Purpose**
            get the standard deviation of the expression data
            (Does not calculate if your data is actually normal)
            
        """
        data = [item["conditions"] for item in self.linearData]       
        data = numpy.array(data)
        return(numpy.std(data))   
