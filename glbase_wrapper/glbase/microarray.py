"""

BUGS
----

TODO
----

* conditions are not always stored as floats... [fixed?]
* Some R functions not implemented in matplotlib yet.
* need to abstract away to draw.
* often the reader will choke on an array... [partially fixed?]
* move the "draw" stuff into a seperate heritable class. [partially implemented]
* to draw, do things like: array.draw.heatmap()
* normaliseToCondition has a bug which is currently bodged over - if the background = 0.0 then it will DivByZero Error.
"""

import sys, os, csv, string, math

from numpy import array
from array import array as qarray
from copy import deepcopy

import config
from flags import *
from genelist import genelist
from draw import draw
if config.MATPLOTLIB_AVAIL:
    import matplotlib.pyplot as plot

class microarray(genelist):
    """
    Container class for microarray data.
    expects array data in a particular format.
    [RefSeq, EntrezGene, Symbol (or Empty), loc (or Empty), condition1, condition2 ...]

    Inherits list-like behaviour from geneList
    """
    def __init__(self, filename, **kargs):
        genelist.__init__(self)

        format = None
        self.name = filename
        self.condition_names = []
        for k in kargs:
            if k == "name":
                self.name = kargs[k]
            if k == "format":
                format = kargs[k] # override the default format.

        if not format:
            format = {"refseq": 2, "entrez": 1, "name": 3,
                    "conditions": {"code": "column[4:]"}, "array_systematic_name": 0,
                    "duplicates_key": False} # use the sniffer first time around?

        self.loadCSV(filename=filename, format=format) # no need for error checking here - it's in genelist now.

        # reopen the file to get the condition headers.
        oh = open(filename, "rU")
        if format.has_key("dialect"):
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

        self.array = {"file": filename, "name": self.name, "hash": self.__hash__()}
        self.cacheLocation = None

        # convert the conditions to floats.
        for item in self:
            conditions = item["conditions"]
            for i, c in enumerate(conditions):
                conditions[i] = float(c)

        self._optimiseData()

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
        Add optional microarray optimisations
        """
        genelist._optimiseData(self) # do the parent

        # generate a serialised version of the array conditions.
        con_names = self.condition_names
        data = {}
        for array_set in self["conditions"]:
            for index, name in enumerate(con_names):
                if not data.has_key(name):
                    data[name] = []
                data[name].append(float(array_set[index])) # get the particular column
        self.serialisedArrayDataDict = data

        # list;
        condition_names = self.getConditionNames()
        self.serialisedArrayDataList = [self.serialisedArrayDataDict[key] for key in condition_names]
        return(True)

    def saveCSV(self, path, filename):
        """
        (Override)
        Some special cases for microarray data.
        save the geneList as a csv
        any geneList saved using saveCSV() you can load it back in using
        loadCSV and the sniffer format - it will always be compatible.
        """
        oh = open(os.path.join(path, filename), "w")
        writer = csv.writer(oh, dialect=csv.excel_tab)

        #format = {"refseq": 0, "entrez": 1, "symbol": 2, "coords": 3,
        #        "conditions": {"code": "column[4:]"}, "array_systematic_name": 1, "duplicates_key": False,
        #        "dialect": csv.excel_tab}

        writeOrder = ["array_systematic_name", "entrez", "refseq", "name"]

        writer.writerow(writeOrder + self.getConditionNames())

        for data in self.linearData:
            line = []
            for key in writeOrder:
                line.append(data[key])
            line = line + data["conditions"] # conditions go last.
            writer.writerow(line)
        oh.close()

    # ----------- overrides/extensions ---------------------------------

    def getArrayName(self): # better as __str__?
        return(self.array["name"])

    def getGenomeName(self):
        return(self.genome.getName())

    def getConditions(self, key, value):
        """
        get the conditions for "key"

        #? this is named wrong?
        It should be something like, get conditions by key?
        I think this is not really what I want to do.
        """
        return(self._findDataByKeyGreedy(key, value)["conditions"])

    def getDataForCondition(self, condition_name):
        """
        get all of the microarray data for a particular condition
        name, returns a list of all the values.
        The list remains in the same order as the overall list,
        so to get e.g. the refseq column, do: l = microarray["refseq"]
        """
        names = self.getConditionNames()
        if not condition_name in names:
            print "Error: No condition named: %s" % condition_name
            return(False)

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

        **Result**

        saves an image to the 'filename' location and
        returns the 'actual filename' that really gets saved (glbase
        will modify the .png ending to .svg or .eps depending upon
        the current setings).
        """
        # checks for option here please.
        assert filename, "you must specify a filename"

        actual_filename = self.draw._heatmap(data=self.serialisedArrayDataDict,
            row_names=self["name"], col_names=self.getConditionNames(),
            filename=filename, **kargs)

        print "Info: Saved the heatmap image: %s" % actual_filename
        return(actual_filename)

    def normaliseToCondition(self, condition_name, bUseFoldChange=True, keep_normed=False):
        """
        normalise all other conditions to condition_name and delete condition name from the array list
        keep_normed will keep the original data, but set it to 0

        This puts the data in the range 0 .. 1 where 0.5 = the same as the control.
        """
        names = self.getConditionNames()

        if condition_name not in names:
            print "Error: condition name: %s is not on this array" % condition_name
            return(None)

        name_index = names.index(condition_name)
        #print name_index

        newl = self.__copy__()
        newl.linearData = []
        newl.condition_names = []

        for item in self:
            #print item
            old_array_data = item["conditions"]
            #print old_array_data
            new_array_data = []
            toNormal = old_array_data[name_index]+0.0000001 # stop divByZero errors.
            for index, datum in enumerate(old_array_data):
                name = names[index] # verbose for clarity
                if name != condition_name:
                    if bUseFoldChange:
                        new_array_data.append(math.pow(2, -(float(toNormal) / float(datum))))
                    else:
                        new_array_data.append(float(toNormal) / float(datum))
                elif keep_normed:
                    new_array_data.append(1.0)

            data_copy = deepcopy(item)
            newl.linearData.append(data_copy)
            data_copy["conditions"] = new_array_data # load the new_array_data over the old one.

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
        print "Info: Normalised to condition: %s" % condition_name
        self._history.append("Normalised To Condition: %s" % condition_name)
        return(newl)

    def drawDotplot(self, x_condition_name, y_condition_name, **kargs):
        """
        draw an X/Y dot plot, get R^2 etc.

        x_condition_name = the name of the er... X condition
        y_condition_name = the name of the er... Y condition

        available key-word arguments:
        xlimits = (0, 2)
        ylimits = (0,2)
        filename = "dotplot.png"
        log = False|True
        """
        # do lots of testing for the existance or not of the condition key
        x_data = self.getDataForCondition(x_condition_name)
        y_data = self.getDataForCondition(y_condition_name)

        if len(x_data) < 100:
            # prefer matplotlib for <100
            if config.MATPLOTLIB_AVAIL:
                usePlotter = "matplot"
            else:
                print "Error: matplotlib not available, cannot drawdotplot()"
                sys.quit()
        else:
            # Use a mapped image instead of this one of r"s smoothScatter
            if config.MATPLOTLIB_AVAIL:
                usePlotter = None
                print "Warning: dotplots for samples > 100 not implemented"

        # defaults:
        xlims = False
        ylims = False
        r_xlims = r.c(0, max(x_data + y_data))
        r_ylims = r.c(0, max(x_data + y_data))
        filename="dotplot.png"
        log = False
        for key in kargs:
            if key == "xlimits": # requires a tuple.
                r_xlims = r.c(kargs[key][0],kargs[key][1])
                xlims = (kargs[key])
            if key == "ylimits": # requires a tuple.
                r_ylims = r.c(kargs[key][0],kargs[key][1])
                ylims = (kargs[key])
            if key == "log":
                if kargs["log"]:
                    pass
                    # log the data.
            if key == "filename":
                filename = kargs[key]

        if usePlotter == "matplot":
            plot.cla()
            plot.scatter(x_data,y_data)
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
            # for the matplot.lib < 100: I want to label everything.
            names = self["name"]
            for i, n in enumerate(names):
                plot.annotate(n, (x_data[i], y_data[i]), size=6)

            plot.savefig(filename)
        print "Info: Saved the dotplot image"
        return(True)

    def drawBoxplot(self, **kargs):
        """
        draw's a boxplot of all conditions.
        """
        ylimits = False
        filename = "Boxplot.png"

        data = self.serialisedArrayDataList
        plot.cla()
        for key in kargs:
            if key == "log":
                if kargs["log"]:
                    data = copy(self.serialisedArrayDataList)
                    for set in data:
                        print set
                        for index, item in enumerate(set):
                            print index, item
                            set[index] = math.log(2, item)
            if key == "ylimits": # requires a tuple.
                ylimits = kargs[key]
            if key == "filename":
                filename = kargs[key]
        # do plot
        plot.boxplot(self.serialisedArrayDataList)
        if ylimits:
            plot.ylim(ylimits)

        plot.xticks(arange(1,len(self.getConditionNames())+1),self.getConditionNames())
        plot.savefig(filename)
        print "Info: Saved the boxplot image: %s" % filename
        return(True)

    def drawDistributionCurves(self, **kargs):
        """
        draws the distributions
        valid keyword arguments:
        filename = filename of the resulting
        window = size of window for moving average
        modifier = undocumented fudge for float based arrays.
        xlimits = a tuple or list of the form: (minimum_x, maximum_y)
        """
        filename = "distribution_plot.png"
        binner_modifier = 100
        window_size = 20
        simple_args = ["filename", "window", "modifier"]
        xlimits = None
        for k in kargs:
            #if k in simple_args:
            #    eval(
            if k == "filename":
                filename = kargs[k]
            if k == "window":
                window_size = kargs[k]
            if k == "modifier": # undocumented fudge for flat arrays.
                binner_modifier = kargs[k]
            if k == "xlimits": # requires a tuple.
                xlimits = kargs[k]

        # normalise data, for each condition.
        plot.cla()
        conditions = self.getConditionNames()
        for c in conditions:
            data = self.getDataForCondition(c)
            a = qarray("L", [0 for x in xrange(int((max(data) + 1)  * binner_modifier))]) # assumes data is 0 -> max bound...
            for d in data:
                a[int(d * binner_modifier)] += 1
            x, n = utils.movingAverage(a, window_size)
            plot.plot(x, n, label=c)
        if xlimits: plot.xlim(xlimits)
        plot.legend()
        plot.savefig(filename)
        print "Info: Saved the Distribustion Curves: %s" % filename
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
            print "Error: Criteria function unavailable."
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

    def _insertCondition(self, condition_name, condition_data, range_bind=None):
        """
        candidate for reveal?
        """
        self.condition_names.append(condition_name)
        max_data = max(condition_data)
        min_data = min(condition_data)
        if len(condition_data) != len(self):
            print "Error: Insertion of array data, wrongly sized"
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

    def sort(self, key):
        """
        This is slightly different from the vanilla genelist's sort - you can pass it the name of
        a condition. Take care to make sure the condition name is not a valid list key.
        The algorithm searches the genelist before searching the array for your particular condition.

        **Arguments**

        key

            must be a valid key in the genelist or the name of an array condition.

        **Result**

        returns True if succesful.

        returns False if no valid.
        """
        if key in self:
            genelist.sort(key) # use the parents sort.
            return(True)
        else:
            names = self.getConditionNames()
            if key in names:
                name_index = names.index(key)
                self.linearData = sorted(self.linearData, cmp=lambda x, y: cmp(x["conditions"][name_index],y["conditions"][name_index])) # the original sort() was overridden.
                self._optimiseData()
                return(True)
        return(False)

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        }
    # define gui interfaces
    #annotate.__gui__ = {"required": {"list": "genelist", "key": "list key", "distance": "int"},
    #    "optional": {"resolution": "filename"}} # bind gui descriptor

    #load.__gui__ = {"required": {"filename": "str"},
    #    "optional": {"format": "option"}} # bind gui descriptor

