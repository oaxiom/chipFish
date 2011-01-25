"""
**Purpose**
    An all-purpose container for microarray data.

**to do**

* conditions are not always stored as floats... [fixed?]
* Some R functions not implemented in matplotlib yet.
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

from numpy import array, arange, meshgrid, zeros, linspace, mean, object_, std
from array import array as qarray # this should be deprecated later.

import config
from flags import *
from genelist import genelist
from draw import draw
from progress import progressbar
from errors import AssertionError, ArgumentError

import matplotlib.pyplot as plot
from pylab import bivariate_normal, griddata # comes from where?
import matplotlib.cm as cm

class microarray(genelist):
    def __init__(self, filename=None, format=None, **kargs):
        """
        **Purpose**

            A container class for microarray data.
            Useful functions for manipulating microarray data in
            relation to other genelists.
            Inherits all genelist methods and implements several microarray
            specific methods.

            microarray analysis in glbase is very underpowered and requires normalised
            and easily manipulatable data. Examples for normalisation are
            genespring output, or output from R and PMA.

        **Arguments**

            filename (Required)

                the filename of the microarray to load.

            format (Optional)

                a format specifier.
                Must include some form of {"conditions": {"code": "column[4:]"}}
                to specifiy the location of the numeric array data.
                glbase will try to guess but will likely fail.

                There are a few special import methods available.
                To import illumina microarray data output from bead studio
                use the special format specifier: "illumina"

        **Returns**

            A microarray instance.
        """
        valig_args = ["filename", "format"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.__init__, k)

        genelist.__init__(self)
        self.name = filename
        self.condition_names = []

        assert filename, "no filename specified"
        assert os.path.exists(os.path.realpath(filename)), "'%s' not found" % filename
        if format:
            assert "conditions" in format, "no 'conditions' specification to collect the data"

        for k in kargs:
            if k == "name":
                self.name = kargs[k]
            if k == "format":
                format = kargs[k] # override the default format.
                assert "conditions" in kargs["format"], "you must provide a 'conditions' entry in the format specifier"

        if not format:
            # this is the export fomrat from GeneSpring 7.
            format = {"refseq": 2, "entrez": 1, "name": 3,
                    "conditions": {"code": "column[4:]"}, "array_systematic_name": 0,
                    "duplicates_key": False} # use the sniffer first time around?

        self.loadCSV(filename=filename, format=format) # no need for error checking here - it's in genelist now.

        # reopen the file to get the condition headers.
        oh = open(filename, "rU")
        if "dialect" in format:
            reader = csv.reader(oh, dialect=format["dialect"])
        else:
            reader = csv.reader(oh)

        self.condition_names = []
        for column in reader:
            exec "names = %s" % format["conditions"]["code"] # yay, more nice happy arbitrary code execution.

            for key in names:
                self.condition_names.append(key)
            break
        oh.close()

        self.conditions = self.condition_names # also accesible
        self.array = {"file": filename, "name": self.name, "hash": self.__hash__()}
        self.cacheLocation = None

        # convert the conditions to floats.
        for item in self:
            conditions = item["conditions"]
            for i, c in enumerate(conditions):
                conditions[i] = float(c)

        self._optimiseData()
        config.log.info("Loaded Microarray data %s items long" % len(self))

    def __repr__(self):
        return("glbase.microarray")

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
        return(self.condition_names)

    def _optimiseData(self):
        """
        (Override)
        (Internal)
        Add microarray optimisations
        """
        genelist._optimiseData(self) # do the parent optimise.

        # generate a serialised version of the array conditions.

        data = {}
        for array_data in self["conditions"]:
            for index, name in enumerate(self.condition_names):
                if not name in data:
                    data[name] = []
                data[name].append(float(array_data[index])) # get the particular column
        self.serialisedArrayDataDict = data
        # this is broken?
        # access aray data by name:
        #self.serialisedArrayDataDict["condition_name"]

        # The __array_data is immutable. serialisedArrayDataDict has been deleted.
        # Instead, in the main list ["array_data"] key points to the row of the numpy array instead.
        #self.__array_data = zeros([len(self["conditions"]), len(self)], dtype=object_) # a numpy version.
        #for r, array_data in enumerate(self["conditions"]):
        #    for c, v in enumerate(array_data):
        #        self.__array_data[r, c] = {"value": v, "stderr": 0.0} # this is the array data type.

        # list;
        self.serialisedArrayDataList = [self.serialisedArrayDataDict[key] for key in self.condition_names]
        # this is just a list version.
        #self.serialisedArrayDataDict[0]
        return(True)

    def saveTSV(self, filename=None, tsv=True, **kargs):
        """
        (Override)
        **Purpose**

            Save the microarray data as a tsv file

        **Arguments**

            filename
                The filename (with a valid path) to save the file to.
                
            tsv (True|False, default=True)

        **Returns**

            returns None
            A saved csv file in filename.
        """
        valig_args = ["filename", "tsv"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.saveCSV, k)

        assert filename, "you must specify a filename"

        oh = open(os.path.realpath(filename), "w")
        if tsv:
            writer = csv.writer(oh, dialect=csv.excel_tab)
        else:
            writer = csv.writer(oh)

        writeOrder = [key for key in self.linearData[0] if key != "conditions"]

        title_row = []
        for k in writeOrder:
            if k in self.linearData[0]: # sample to see if we have this key
                title_row.append(k) # should mimic below
        writer.writerow(title_row + self.getConditionNames())

        for data in self.linearData:
            line = []
            for key in writeOrder:
                if key in data:
                    line.append(data[key])
            writer.writerow(line + data["conditions"])# conditions go last.
        oh.close()
        config.log.info("Saved a csv file to '%s'" % filename)
        return(None)

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

            return a copy of the microarray, but only containing
            the condition names specified in conditions

        **Arguments**

            conditions (Required)

                A list, or other iterable of condition names to extract
                from the microarray. Every condition name must be present
                on the microarray.

        **Result**

            A new microarray object with the same settings as the original,
            but containing only the microarray conditions specified in
            the 'conditions' argument.
        """
        valig_args = ["conditions"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.sliceConditions, k)

        assert conditions, "You must specify a list of conditions to keep"
        assert isinstance(conditions, list), "You must specify a list of conditions to keep"
        for item in conditions:
            assert item in self.condition_names, "'%s' condition not found on this microarray" % item

        newl = copy.deepcopy(self)

        copymask = [] # work out a mask to extract the correct array columns.
        for name in conditions:
            for i, c in enumerate(self.condition_names):
                if c == name:
                    copymask.append(i)

        for index, item in enumerate(self.linearData):
            new_array_data = []
            for c in copymask: # use the mask to collect the array entries.
                new_array_data.append(item["conditions"][c])
            newl.linearData[index]["conditions"] = new_array_data

        newl._history.append("sliced conditions, kept %s" % "".join(conditions))
        newl.condition_names = conditions
        newl._optimiseData()
        return(newl)

    def getDataForCondition(self, condition_name):
        """
        get all of the microarray data for a particular condition
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
        for item in con:
            l.append(item[index])
            #print con[index]
        return(l)

    def drawHeatmap(self, filename=None, **kargs):
        """
        **Purpose**

        draw a simple heatmap of the current microarray data.

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
        for index, name in enumerate(self.condition_names):
            newdata[name] = data[index] # get the particular column

        actual_filename = self.draw._heatmap(data=newdata,
            row_names=self["name"], col_names=self.getConditionNames(),
            filename=filename, **kargs)

        config.log.info("Saved the heatmap image: %s" % actual_filename)
        return(actual_filename)

    def normaliseToCondition(self, condition_name, bUseFoldChange=True, keep_normed=False, **kargs):
        """
        **Purpose**
            normalise all other conditions to condition_name and delete
            condition name from the array list

        **Arguments**
            keep_normed (boolean, Optional)
                keep the original data, but set it to 0

            ??? not true right:
            This puts the data in the range 0 .. 1 where 0.5 = the same as the control.

        **Returns**
            returns the newly normalised list
        """
        names = self.getConditionNames()

        assert condition_name in names, "condition name: %s is not on this array" % condition_name

        name_index = names.index(condition_name)
        #print name_index

        newl = self.__copy__()
        newl.linearData = []
        newl.condition_names = []

        p = progressbar(len(self.linearData))
        for index, item in enumerate(self.linearData):
            #print item
            old_array_data = item["conditions"]
            #print old_array_data
            new_array_data = []
            toNormal = old_array_data[name_index]+0.0000001 # stop divByZero errors.
            for i, datum in enumerate(old_array_data):
                name = names[i] # verbose for clarity
                if name != condition_name:
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
            newl.condition_names = names
        else:
            # delete the old label
            newl.condition_names = []
            for name in names:
                if name != condition_name:
                    newl.condition_names.append(name)
                elif keep_normed:
                    newl.condition_names.append(name)

        newl._optimiseData()
        config.log.info("Normalised to condition: %s" % condition_name)
        self._history.append("Normalised To Condition: %s" % condition_name)
        return(newl)

    def drawDotplot(self, x_condition_name, y_condition_name, filename=None, **kargs):
        """
        **Purpose**
            draw an X/Y dot plot, get R^2 etc.

        **Arguments**
            x_condition_name = the name of the er... X condition
            y_condition_name = the name of the er... Y condition

            available key-word arguments:
            xlimits = (0, 2)
            ylimits = (0,2)
            filename = "dotplot.png"
            log = False|True

        **Returns**
            the actual filename saved as and a new image in filename.
        """
        assert filename, "no filename specified"
        assert x_condition_name in self.serialisedArrayDataDict, "%s x-axis condition not found" % x_condition_name
        assert y_condition_name in self.serialisedArrayDataDict, "%s y-axis condition not found" % y_condition_name

        x_data = self.getDataForCondition(x_condition_name)
        y_data = self.getDataForCondition(y_condition_name)

        if "log" in kargs and kargs["log"] == "log2":
            x_data, y_data = self.__log_transform_data([x_data, y_data], 2)

        # defaults:
        xlims = [min(x_data + y_data), max(x_data + y_data)]
        ylims = [min(x_data + y_data), max(x_data + y_data)]
        log = False

        for key in kargs:
            if key == "log":
                if kargs["log"]:
                    pass
                    # log the data.

        plot.cla()
        plot.subplot(111)

        plot.scatter(x_data, y_data, c="blue", alpha=0.3, edgecolors='none')

        plot.xlabel(x_condition_name)
        plot.ylabel(y_condition_name)
        if xlims:
            plot.xlim(xlims)
        else:
            plot.xlim((min(x_data), max(x_data + y_data)))
        if ylims:
            plot.ylim(ylims)
        else:
            plot.ylim((min(y_data), max(y_data + x_data)))

        if len(x_data) < 100:
            # for the matplot.lib < 100: I want to label everything.
            names = self["name"]
            for i, n in enumerate(names):
                plot.annotate(n, (x_data[i], y_data[i]), size=6)

        plot.savefig(filename)
        real_filename = filename # not supported yet.
        config.log.info("Saved the dotplot image '%s'" % real_filename)
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
        actual_filename = self.draw._boxplot(data=data, filename=filename, xticklabels=self.getConditionNames())

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
                data.append(row)
            return(data)
        else:
            return(serialisedArrayDataList)

    def drawCurves(self, filename=None, **kargs):
        """
        **Purpose**

        draw a bell-curve diagram of the array expression.

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
                
            verbose (True|Fals, default=False
                print out the means and standard distributions.

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
        plot.cla()

        if "verbose" in kargs and kargs["verbose"]:
            print "name\tmean\tstd" 

        for i, c in enumerate(self.getConditionNames()):
            plot.hist(data[i], bins=window_size, histtype="step", label=c, **extra_args)
            m =  mean(data[i])
            d = std(data[i])
            plot.axvline(x=m, color="red")
            plot.axvline(x=m-d, color='grey', ls=":")
            plot.axvline(x=m+d, color='grey', ls=":")
            if "verbose" in kargs and kargs["verbose"]:
                print "%s\t%.2f\t%.2f" % (c, m, d)

        if xlimits: plot.xlim(xlimits)
        plot.legend()
        plot.savefig(filename)
        config.log.info("Saved the Distribustion Curves: %s" % filename)
        return(True)

    def getDataByCriteria(self, **kargs):
        """
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
        self.condition_names.append(condition_name)
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
            A new microarray object with probes that fail to pass removed.
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

    def sort(self, key):
        """
        This is slightly different from the vanilla genelist's sort - you can pass it the name of
        a condition. Take care to make sure the condition name is not a valid list key.
        The algorithm searches the genelist before searching the array for your particular condition.

        Also take care with this one: It is one of the few in-place list
        modifiers.

        **Arguments**

        key

            must be a valid key in the genelist or the name of an array condition.

        **Result**

        returns True if succesful.

        returns False if no valid.
        """
        assert (key in self.linearData[0]) or key in self.getConditionNames(), "'%s' search key not found in list or array data" % key

        if key in self.linearData[0]:
            return(genelist.sort(self, key)) # use the parents sort.
        else:
            names = self.getConditionNames()
            if key in names:
                name_index = names.index(key)
                self.linearData = sorted(self.linearData, cmp=lambda x, y: cmp(x["conditions"][name_index],y["conditions"][name_index])) # the original sort() was overridden.
                self._optimiseData()
                return(True)
        return(False)

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
        assert key in self, "key '%s' not found in the microarray" % key

        mapped_scores = []

        plot.cla()
        fig = plot.figure()
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










