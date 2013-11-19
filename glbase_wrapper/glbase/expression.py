"""
**Purpose**
    An all-purpose container for transcript expression data.

**to do**

* often the reader will choke on an array... [partially fixed?]
* scipy.cluster relies of recursion to perform the clustering.
    On large sized arrays it chokes horribly.
"""

import sys, os, csv, string, math, copy

from operator import itemgetter

import numpy, scipy
from numpy import array, arange, meshgrid, zeros, linspace, mean, object_, std
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from scipy.cluster.vq import vq, kmeans, whiten, kmeans2
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plot
from pylab import bivariate_normal, griddata # comes from where?
import matplotlib.cm as cm
from scipy.stats import ttest_ind, pearsonr
from scipy import stats

import config, utils
from flags import *
from base_expression import base_expression
from draw import draw
from progress import progressbar
from errors import AssertionError, ArgumentError
from genelist import genelist
from location import location
from pca import pca

class expression(base_expression):
    def __init__(self, loadable_list=None, filename=None, format=None, expn=None, **kargs):
        """
        **Purpose**
            The base container for expression data.

            Please not though that:
            Expression analysis in glbase requires normalised
            and easily manipulatable data. Examples for normalisation are
            genespring output, or output from R and PMA, Cufflinks RSEM etc...

        **Arguments**

            filename (Required, one of loadable_list or filename)

                the filename of the microarray to load.

            loadable_list (Required, one of loadable_list or filename)
            
                a genelist-like object I can use to construct the list from. If you use this then
                expn should be a list of keys to extract the expression data from.
                
            expn (Required)
                If filename:
            
                Some sort of descriptor telling me where the expression data actually is.
                For example:
                    "column[3:]" (Take each column from column3 onwards, until the end column)
                Or:
                    "column[4::2]" (Take every second item, from column 4 until the end)
                Or:
                    "column[5:8]" (Take column 5 through 7 - 7 is not a typo. The lists are 
                        zero-ordered and closed)
                        
                In fact you can even give compound statements:
                    "[column[7], column[17]]" (Take columns 7 and 17)
                    
                The only rules are it must be a valid piece of python code.
                
                If a loadable_list:
                
                This should be the name of the keys used to extract the expresion data

            err (Optional)
                Some sort of descriptor for where to get the error data from.
                NOT IMPLEMENTED
            
            cv_err (Optional)
                Some sort of descriptor for where to get confidence intervals from.
                This should be a tuple, as confidence intervals are 'lo' and 'hi'.
                NOT IMPLEMENTED

            cond_names (Optional)
                A list of the condition names (in order) if glbase is not working them 
                out itself.
            
            name (Optional)
                By default expression will use the filename (removing any .txt, .tsv, etc) from
                the ends of the names      
                
            silent (Optional, default=False)
                Do not output any reports (primarily this is for internal functions)
        """
        if loadable_list:
            base_expression.__init__(self, loadable_list=loadable_list, expn=expn, **kargs)
        else:
            base_expression.__init__(self, filename=filename, expn=expn, format=format, **kargs)
            
        self.pca = None # Store for PCA object

    def __repr__(self):
        return("glbase.expression")

    def __getitem__(self, index):
        """
        Confers:
        
        a = expn["condition_name"]
        
        and inherits normal genelist slicing behaviour
        """
        if index in self._conditions:
            return(self.getDataForCondition(index))
        return(base_expression.__getitem__(self, index)) # otherwise inherit

    def sort_sum_expression(self):
        """
        sort by the sum of conditions

        **Arguments**
            None
            
        """
        self.linearData = sorted(self.linearData, cmp=lambda x, y: cmp(sum(x["conditions"]), sum(y["conditions"]))) # the original sort() was overridden.
        self._optimiseData()
        return(True)

    def multi_sort(self, keys):
        """
        **Purpose**
            Sort a genelist using multiple keys. 
            
            This version is tailored for expression objects, and will accept the name of a condition
            
        **Arguments**
            keys (Required)
                A list of key names to sort. Sorting will first sort keys[0] then key[1] through key[n]
                
        **Returns**
            returns True if it completes.
            sorts the list IN PLACE.
        """
        #assert key, "No such key '%s'" % key
        #assert key in self.linearData[0], "Data does not have key '%s'" % key
        
        comparers = []
        for k in keys:
            if k in self._conditions:
                kind = self._conditions.index(k)
                comparers.append(lambda x: x["conditions"][kind])
            else:
                comparers.append(itemgetter(k))
        
        def comparer(left, right):
            for fn in comparers:
                result = cmp(fn(left), fn(right))
                if result:
                    return result
            else:
                return 0
        
        self.linearData = sorted(self.linearData, cmp=comparer)
        self._optimiseData()
        return(True)

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

    # ----------- overrides/extensions ---------------------------------

    def getGenomeName(self):
        if not self.genome:
            return("Genome not bound")
        return(self.genome.getName())

    def getColumns(self, return_keys=None, strip_expn=False):
        """
        **Purpose**
            return a new expression only containing the columns specified in return_keys (a list)
            
            This version for expression objects will preserve the conditions and expresion data.
            
        **Arguments**
            return_keys (Required)
                A list of keys to keep
            
            strip_expn (Optional, default=False)
                If True then remove the expression and err keys (if present).
        
        **Returns**
            A new expression object or if strip_expn=True then a genelist
        """
        assert isinstance(return_keys, list), "return_keys must be a list"

        if strip_expn:
            newl = genelist()
            newl.name = str(self.name)
        else:
            newl = self.shallowcopy()
            newl.linearData = []
            if not "conditions" in return_keys and "conditions" in self.linearData[0]:
                return_keys.append("conditions")
            if not "err" in return_keys and "err" in self.linearData[0]:
                return_keys.append("err")           

        for item in self.linearData:
            newd = {} # Could be done with dict comprehension.
            for key in return_keys:
                newd[key] = item[key] # This is wrong? It will give a view? 
                
            newl.linearData.append(newd)
        newl._optimiseData()
        
        config.log.info("getColumns(): got only the columns: %s" % ", ".join(return_keys))
        return(newl)

    def strip_errs(self):
        """
        **Purpose**
            Remove the any err keys if present. 
            Sometimes the err keys can be preserved inappropriately (Mostly due to running it through
            a method which does not support rescaling errors). This function removes the keys from the list
            
            This is an IN PLACE method.
            
        **Arguments**
            None
            
        **Returns**
            None
        """
        for item in self.linearData:
            del item["err"]
        self._optimiseData() # I think this does nothing at the moment, but just in case I ever fix err key handling
        return(None)

    def sliceConditions(self, conditions=None, **kargs):
        """
        **Purpose**

            return a copy of the expression-data, but only containing
            the condition names specified in conditions
            
            Note that you can also use this method to change the order of the conditions.
            Just slice all of the keys in the order you want them to occur in.
            
            Additionally, you can use this to replicate a key, e.g.
            
            gl = gl.sliceConditions(["cond1", "cond2", "cond1"])
            
            will now give you an expression set with two 'cond1' conditions

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
        assert isinstance(conditions, list) or isinstance(conditions, set), "You must specify a list or set of conditions to keep"
        for item in conditions:
            assert item in self._conditions, "'%s' condition not found in this expression-data" % item

        newgl = self.shallowcopy()
        newgl.linearData = []

        copymask = [] # work out a mask to extract the correct array columns.
        for name in conditions:
            for i, c in enumerate(self._conditions):
                if c == name:
                    copymask.append(i)
                    break

        for index, item in enumerate(self.linearData):
            newi = copy.deepcopy(item) # safe copy available
            newi["conditions"] = [copy.deepcopy(item["conditions"][c]) for c in copymask]
            if "err" in item:
                newi["err"] = [copy.deepcopy(item["err"][c]) for c in copymask]                   
            
            newgl.linearData.append(newi)
            
        newgl._history.append("sliced conditions, kept %s" % "".join(conditions))
        newgl._conditions = conditions
        newgl._optimiseData()
        config.log.info("sliceConditions(): sliced for %s conditions" % (len(newgl[0]["conditions"]),))
        return(newgl)

    def getDataForCondition(self, condition_name):
        """
        **Purposse**
            get all of the expression-data data for a particular condition
            name, returns a list of all the values.
            The list remains in the same order as the overall list,
        """
        #print self.serialisedArrayDataDict.keys()
        assert condition_name in self.getConditionNames(), "No condition named '%s' in this expression object" % condition_name

        return(self.serialisedArrayDataDict[condition_name])

    def getExpressionTable(self):
        """
        **Purpose**
            Return the entire expression table. Note that rows and columns are not labelled.
            Will return a numpy list
            
        **Arguments**
            None
        """
        ll = []
        for item in self.linearData:
            ll.append(item["conditions"])
        return(numpy.array(ll))        

    def coerce(self, new_type):
        """
        **Purpose**
            Semi-internal/obscure function. Coerces the data in condition into
            the type specified by new type. Primarily this is to convert the 
            expression data from/to integers or floats for downstream R problems.
            
        **Arguments**
            new_type (Required)
                generally int or float
            
        **Returns**
            None
            THIS IS AN IN-PLACE CONVERSION
        """
        for item in self.linearData:
            item["conditions"] = [new_type(i) for i in item["conditions"]]
        return(None)

    def drawHeatmap(self, filename=None, row_label_key="name", **kargs):
        """
        Obselete alias for heatmap()
        """
        return(self.heatmap(filename=filename, row_label_key=row_label_key, **kargs))
    
    def heatmap(self, filename=None, row_label_key="name", **kargs):
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

            row_tree (Optional, default=False)
                provide your own tree for drawing. Should be a Scipy tree. row_labels and the data
                will be rearranged based on the tree, so don't rearrange the data yourself.

            col_tree (Optional, default=False)
                provide your own tree for ordering the data by. See row_tree for details.
                This one is applied to the columns.

            highlights (Optional, default=None)
                sometimes the row_labels will be suppressed as there is too many labels on the plot. 
                But you still want to highlight a few specific genes/rows on the plot.
                Send a list to highlights that matches entries in the row_names.

            discretize (Optional, default=False)
                change the colourmap (either supplied in cmap or the default) into a 'discretized' version
                that has large blocks of colours, defined by the number you send to discretize

            col_norm (Optional, default=False)
                normalise each column of data between 0 .. max => 0.0 .. 1.0
                
            row_norm (Optional, default=False)
                similar to the defauly output of heatmap.2 in R, rows are normalised 0 .. 1
                                
            row_font_size (Optional, default=guess suitable size)
                the size of the row labels (in points). If set this will also override the hiding of
                labels if there are too many elements. 
                
            col_font_size (Optional, default=8)
                the size of the column labels (in points)
            
            heat_wid (Optional, default=0.25)
                The width of the heatmap panel. The image goes from 0..1 and the left most
                side of the heatmap begins at 0.3 (making the heatmap span from 0.3 -> 0.55).
                You can expand or shrink this value depending wether you want it a bit larger
                or smaller.
                
            heat_hei (Optional, default=0.85)
                The height of the heatmap. Heatmap runs from 0.1 to heat_hei, with a maximum of 0.9 (i.e. a total of 1.0)
                value is a fraction of the entire figure size.
                
            colbar_label (Optional, default="expression")
                the label to place beneath the colour scale bar            
        **Result**

            saves an image to the 'filename' location and
            
            returns a dictionary containing the real_filename and re-ordered labels after clustering (if any).
            
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
        for index, name in enumerate(self._conditions):
            newdata[name] = data[index] # get the particular column

        res = self.draw.heatmap(data=newdata,
            row_names=self[row_label_key], col_names=self.getConditionNames(),
            filename=filename, **kargs)
            
        config.log.info("heatmap(): Saved %s" % res["real_filename"])
        return(res)

    def __fold_change(self, c1v, c2v, log=2):
        """
        Calculate the fold-change of two values
        assumes values have already been padded.
        """
        try:
            if c2v > c1v:
                if log:
                    return(math.log((c2v/c1v), log))
                else:
                    return((c2v/c1v))
            else:
                if log:
                    return(-math.log(c1v/c2v, log))
                else:
                    return(-(c1v/c2v))
                                                
        except (OverflowError, ZeroDivisionError, ValueError):
            if c2v > c1v:
                config.log.error("(%.2f/%.2f) failed" % (c2v, c1v))
            else:
                config.log.error("(%.2f/%.2f) failed" % (c1v, c2v))
            raise Exception("__fold_change() encountered an error, possibly the pad value is too small, or you are trying to apply fold-change to log transformed data")

    def normaliseToCondition(self, condition_name, use_fold_change=False, keep_normed=False, pad=1e-6, 
        log=2, **kargs):
        """
        **Purpose**
            normalise all other conditions to condition_name and delete
            condition name from the expression-data list

        **Arguments**
            keep_normed (boolean, Optional)
                keep the original data, but set it to 0

            use_fold_change (Optional, default=False)
                by default glbase performs:
                    2^ -[(control_sample) / (experimental_sample)] 
                this will result in values bound between 0...1 with no change being 0.5.
                
                If use_fold_change is set to True a value will be calulated such that 
                positive values are the fold-change upwards and -ve values are fold-change down
                wards from the control_sample

            pad (Optional, default=1.0E-06)
                The amount to pad the divisor so as to avoid divide by zero errors and
                overflow errors.
        
            log (Optional, default=2)
                log transform the 'fold-change' by the base defined in log (usually 2, so that logFC=1 = FC=2)
                Only used if use_fold_change = True. Set to False if you don't want to log transform FC.

        **Returns**
            returns the newly normalised list
        """
        names = self.getConditionNames()

        assert condition_name in names, "condition name: %s is not in this expression-data" % condition_name

        name_index = names.index(condition_name)
        #print name_index

        newl = self.shallowcopy()
        newl.linearData = []
        newl._conditions = []

        p = progressbar(len(self.linearData))
        for index, item in enumerate(self.linearData):
            #print item
            old_array_data = item["conditions"]
            #print old_array_data
            new_array_data = []
            toNormal = old_array_data[name_index]+pad # stop divByZero errors.
            for i, datum in enumerate(old_array_data):
                if names[i] != condition_name:
                    if use_fold_change:
                        c1v = toNormal
                        c2v = datum+pad
                        #print c2v, c1v, (c2v/c1v), math.pow(2, -(float(toNormal) / (float(datum)+pad)))
                        new_array_data.append(self.__fold_change(c1v, c2v, log=log))
                    else:
                        new_array_data.append(math.pow(2, -(float(toNormal) / (float(datum)+pad))))
                elif keep_normed:
                    if use_fold_change:
                        new_array_data.append(0.0)
                    else:
                        new_array_data.append(0.5)

            data_copy = copy.deepcopy(item)
            newl.linearData.append(data_copy)
            data_copy["conditions"] = new_array_data # load the new_array_data over the old one.

            p.update(index)

        # rebuild the condition names (optimiseData can't handle this)
        if keep_normed:
            newl._conditions = names
        else:
            # delete the old label
            newl._conditions = []
            for name in names:
                if name != condition_name:
                    newl._conditions.append(name)
                elif keep_normed:
                    newl._conditions.append(name)

        newl._optimiseData()
        config.log.info("normaliseToCondition(): '%s'" % condition_name)
        return(newl)

    def norm_multi_fc(self, conds=None, pad=1e-6, log=2, **kargs):
        """
        **Purpose**
            transform the data, but do it relative to multiple base lines.
            
            For example suppose you had the following samples:
            
            sampleA, sampleB, sampleC, sampleD, sampleE
            
            And sy you wanted to transform the expression values into fold-change, but
            sampleB should be transformed relative to sampleA whilst sampleD and E are relative to sampleC.
            
            You'd do this:
            
            newexpn = expn.norm_multi_fc({"sampleA": ["sampleB"], "sampleC": ["sampleD", "sampleE"]})
            
            mean_replicates is now designed to estimate the error stored in "err" key if available.
            
            The estimate is not guaranteed to be accurate though so use with care. Unfortunately the expression class
            cannot store assymetric errors, so the max() of the error is instead taken. 
            
        **Arguments**
            conds (Required)
                A dictionary, in the form:
                    {"norm_sample": [A, B, list of samples to normalise to norm_sample],
                    "another_norm_sample": [C, D, other samples],
                    ...
                    }
                Each sample in the lists [] will be normalised to the sample named as the key in the dict.
                
            pad (Optional, default=1e-6)
                A pad value to avoid Infinity errors and DivideByZero errors.
                
            log (Optional, default=2)
                Base to log transform the fold-change by. Set to None/False  if you don't want to
                log transform the data.
        **Returns**
            A new expression object.
        """
        assert conds.keys(), "conds must be a dictionary object"
        #if "err" in self.linearData[0]:
        #    config.log.warning("'err' key in this data, norm_multi_fc() will not correct the errors, I will delete this key")
        
        newl = []
        newe = None
        for item in self.linearData:
            newc = [0 for n in item["conditions"]]
            for k in conds:
                kin = self._conditions.index(k)
                newc[kin] = 0.0
                    
                for n in conds[k]:
                    nin = self._conditions.index(n)
                    newc[nin] = self.__fold_change(item["conditions"][kin]+pad, item["conditions"][nin]+pad, log=log)
                    if "err" in newc:
                        del newc["err"]
                        
            if "err" in item:
                newe = [0 for n in item["err"]]
                for k in conds:
                    kin = self._conditions.index(k)
                    newe[kin] = 0.0 # by definition.
                    
                    for n in conds[k]:
                        nin = self._conditions.index(n)
                        expn_lo_err = item["conditions"][nin] - item["err"][nin] + pad
                        expn_base = item["conditions"][kin] + pad # I ignore any potential error here.
                        expn_hi_err = item["conditions"][nin] + item["err"][nin] + pad
                        up = self.__fold_change(expn_base, expn_lo_err, log=log)
                        dn = abs(self.__fold_change(expn_base, expn_hi_err, log=log))
                        newe[nin] = max(up, dn)
                        #print (up,dn)
            
            newi = item.copy()
            newi["conditions"] = newc
            if newe and "err" in newe:
                newi["err"] = newe
            newl.append(newi)
        
        newgl = self.shallowcopy()
        newgl.linearData = newl
        newgl._optimiseData()
        return(newgl)        

    def mean_replicates(self, *reps, **kargs):
        """
        **Purpose**
            replace replicates with a single column representing the mean of the replicates.
            
            Don't send the same condition name as part of two other replicates. this function 
            will fail silently.
            
            The standard error will be stored in a new key, 'err' in the same order as in each 
            condition. This key is used by some downstream drawing tools for error bars.
            
        **Arguments**
            Accepts a single (unnamed) argument and keywords (see below)
                The condition names of replicates to merge, as a set of lists of replicates
                specifying the name of the conditions to merge as replicates. e.g.:
                
                newexpn = expn.mean_replicates(["condA_rp1", "condA_rp2"], ["condB_rp1", "condB_rp2", "condB_rp3"])
                
            Keyword args:
            threshold (Optional, default=0.8)
                by default mean_replicates will measure the Pearson correlation between samples and 
                if the score is below threshold will emit a warning
                
            output_pears (Optional, default=False)
                If a filename, mean_replicates will output a square table for all of the Pearson
                correlation values for each pair of RNA samples, where replicates are available.
        
        **Returns**
            A new expression object containing mean-ed values for each set of replicates.
        """
        
        newe = []
        all_conds = self._conditions
        all_reps = set([x for sublist in reps for x in sublist]) # I still don't know how this works.
        done = []
        
        threshold = 0.8
        if "threshold" in kargs and kargs["threshold"]:
            threshold = kargs["threshold"]
        output_pears = False
        if "output_pears" in kargs and kargs["output_pears"]:
            output_pears = kargs["output_pears"]
            
        for index, item in enumerate(self.linearData):
            newi = copy.deepcopy(item)
            newi["conditions"] = []
            newi["err"] = [] # blank any preexisting keys
            done = []
            new_cond_names = []
            for cind, cond in enumerate(self._conditions): # preserve order of conditions
                if cond in all_reps: # Its one of the reps to merge
                    # check it's not done already:
                    if cond not in done:
                        # get the p it is in:
                        p = [i for i in reps if cond in i][0] # the names of the set of replicates to merge
                        pinds = [self._conditions.index(i) for i in p]
                        expn_vals = [item["conditions"][i] for i in pinds]
                        mean = sum(expn_vals) / float(len(pinds))
                        err = numpy.std(expn_vals) / math.sqrt(len(expn_vals))
                        newi["conditions"].append(mean)
                        newi["err"].append(err)
                        
                        # add all reps to the done list:
                        [done.append(i) for i in p]
                        new_cond_names.append(p[0]) # merge into the 0th key
                else: # not to be merged or modified, so just add it to conditions. 
                    newi["conditions"].append(item["conditions"][cind]) 
                    newi["err"].append(0.0)
                    new_cond_names.append(cond) 
            newe.append(newi) # put it on the new expn object

        pear_out = numpy.zeros([len(self._conditions), len(self._conditions)])         
        # check pairs for pearson correlation 
        for r in reps:
            # r can be a list n entries long. I need to compare all vs all
            for i1, p1 in enumerate(r):
                for i2, p2 in enumerate(r):
                    if i1 != i2 and i1 < i2:
                        p1d = self.getDataForCondition(p1)
                        p2d = self.getDataForCondition(p2)
                        corr = scipy.stats.pearsonr(p1d, p2d)
                        if corr[0] < threshold:
                            config.log.warning("Samples '%s' vs '%s', pearson=%.2f" % (p1, p2, corr[0]))
                        if output_pears:
                            p1ind = self._conditions.index(p1)
                            p2ind = self._conditions.index(p2)
                            pear_out[p1ind, p2ind] = corr[0]
                            pear_out[p2ind, p1ind] = corr[0]
                    elif i1 == i2 and output_pears:
                        p1ind = self._conditions.index(p1)
                        pear_out[p1ind, p1ind] = 1.0

        if output_pears:
            oh = open(output_pears, "w")
            oh.write("\t%s\n" % "\t".join(self._conditions))
            for i, row in enumerate(pear_out):
                d_out = []
                for d in row:
                    if d == 0:
                        d_out.append("")
                    else:
                        d_out.append(str(d))
                oh.write("%s\t%s\n" % (self._conditions[i], "\t".join(d_out)))
            oh.close()
            config.log.info("mean_replicates(): Saved Pearson correlation table '%s'" % (output_pears,))
                        
        newgl = self.shallowcopy()
        newgl.linearData = newe
        newgl._conditions = new_cond_names
        newgl._optimiseData()
        config.log.info("mean_replicates(): Started with %s conditions, ended with %s" % (len(self[0]["conditions"]), len(newgl[0]["conditions"])))
        return(newgl)               

    def add_fc_key(self, key="fc", cond1=None, cond2=None, log=2, pad=1.0E-06, and_err=False, **kargs):
        """
        **Purpose**
            Add in a fold-change key for the fold change from cond1 to cond2.
            
            Note that it will pad the values by 1.0E-08 to avoid divide by zero errors.
            
            You can change this value with the 'pad' argument if it is too large/small
            
        **Arguments**
            key (Optional, default="fc")
                The key name to load the fold-change value into
                
            cond1 (Required)
                condition name 1
                
            cond2 (Required)
                condition name 2
                
            pad (Optional, default=1.0E-06)
                The amount to pad the divisor so as to avoid divide by zero errors and
                overflow errors.
                
            log (Optional, default=2)
                By default fold-changes are log2 transformed. This means a fold-change of 
                2.0 becomes 1.0 and no change is no 0.0. Set this to None to disable this 
                behaviour
            
            and_err (Optional, default=False)
                and estimate the err and load into a key name as specified by and_err.
                Note that your expression object must have an 'err' key to estimate the error.
                Also, expresion objects can't hold assymetric error bars. Hence the minimum value
                will be taken. Now, that may sound a bit odd, but as fold-changes are commonly
                log transformed taking the maximum value will overestimate the error, whilst the min 
                value will show an accurate 'upward' error bar and a somewhat inaccurate downward error
                bar.
                 
        **Returns**
            A new expression object containing the <fc> key.
        """
        assert cond1, "add_fc_key(): you must specify a name for condition1"
        assert cond1 in self._conditions, "add_fc_key(): appears '%s' not in this expression-object" % cond1
        assert cond2, "add_fc_key(): you must specify a name for condition1"
        assert cond2 in self._conditions, "add_fc_key(): appears '%s' not in this expression-object" % cond2
        if and_err:
            assert "err" in self.linearData[0], "add_fc_key(): 'err' key not found in this genelist"
        
        newl = self.deepcopy() # full copy
        
        c1i = self._conditions.index(cond1)
        c2i = self._conditions.index(cond2)
        
        for item in newl:
            c1v = item["conditions"][c1i]+pad
            c2v = item["conditions"][c2i]+pad
            try:
                item[key] = self.__fold_change(item["conditions"][c1i]+pad, item["conditions"][c2i]+pad, log=log)
            except OverflowError, DivByZeroError:
                if c2v > c1v:
                    config.log.error("(%.2f/%.2f) failed" % (c2v, c1v))
                else:
                    config.log.error("(%.2f/%.2f) failed" % (c1v, c2v))
                raise Exception("add_fc_key(): encountered an error, possibly the pad value '%s' is too small, or you need to log transform the data as you are getting an Infinity result" % pad)
            
            if and_err: # If it got here then the above try probably passed.
                expn_lo_err = item["conditions"][c2i] - item["err"][c2i] + pad
                expn_base = item["conditions"][c1i] + pad # I ignore any potential error here.
                expn_hi_err = item["conditions"][c2i] + item["err"][c2i] + pad
                up = abs(item[key] - (self.__fold_change(expn_base, expn_lo_err, log=log)))
                dn = abs(item[key] - self.__fold_change(expn_base, expn_hi_err, log=log))
                item[and_err] = min(up,dn)
                #print item["name"]
                #print item["conditions"]
                #print item["err"]
                #print expn_lo_err, expn_base, expn_hi_err, up, dn, item[key]
            
        newl._optimiseData()
        config.log.info("add_fc_key(): Added fold-change key '%s'" % key)
        return(newl)

    def filter_by_fc(self, fckey=None, direction="any", value=2.0, **kargs):
        """
        **Purpose**
            Filter data which passes some sort of fold-change, defined in a previous key generated
            by add_fc_key()
            
        **Arguments**
            fckey (Required)
                The name of the fc_key to use, previously generated by add_fc_key()
                
            direction (Optional, default="any", values=["up", "down", "dn", "any"])
                The direction of change. "up" is +value, "down"/"dn" is -value and "any" is
                either +value or -value.
                
            value (Optional, default=2.0)
                The value of change required to pass the test.
                comparisons are evaluated as >= (i.e. greater than or equal to)  
                
        **Returns**
            A new expression-object containing only the items that pass.    
        """
        assert fckey, "filter_by_fc(): 'fckey' argument is required"
        assert fckey in self.linearData[0], "filter_by_fc(): '%s' not found in this expression object" % fckey
        assert direction in ("up", "down", "dn", "any"), "filter_by_fc(): direction argument '%s' not recognised" % direction
        
        new_expn = []
        
        for item in self.linearData:
            if direction == "up":
                if item[fckey] >= value:
                    new_expn.append(item)
                    
            elif direction in ("down", "dn"):
                if item[fckey] <= -value:
                    new_expn.append(item)
                    
            elif direction == "any":
                if item[fckey] >= value or item[fckey] <= -value:
                    new_expn.append(item)
                
        ret = expression(loadable_list=new_expn, cond_names=self._conditions)
        
        rep_d = {"up": "+", "dn": "-", "down": "-", "any": "+\-"}
        
        config.log.info("filter_by_fc(): Filtered expression by fold-change '%s' %s%s, found: %s" % (fckey, rep_d[direction], value, len(ret)))
        return(ret)         

    def filter_by_value(self, key=None, evaluator=None, value=None, **kargs):
        """
        **Purpose**
            Filter data based on a key with some numeric data.
            
            If you skip the keyword arguments you can write nice things like this:
            
            newdata = expn.filter_by_value("q-value", "<", 0.05)
            
        **Arguments**
            key (Required)
                The name of the key to use to filter the data on.
                or a name of a condition in the expression object to filter on.
                
            evaluator (Required, values=["gt", "lt", "gte", "lte", "equal"])
                The comparator. 
                    gt = '>' greater than value
                    lt = '<' less than value
                    gte = '>=' greater than or equal to value
                    lte = '<=' less than or equal to value
                    equal = "==" equal to value
                    
                    You can also send the > < >= <= or == as a string symbol as well.
                                
            value (Required)
                The value of change required to pass the test.
                
        **Returns**
            A new expression-object containing only the items that pass.    
        """
        assert key, "filter_by_value(): 'key' argument is required"
        assert evaluator in ("gt", "lt", "gte", "lte", "equal", ">", "<", ">=", "<=", "=="), "filter_by_value(): evaluator argument '%s' not recognised" % evaluator
        
        if key in self._conditions:
            its_a_condition = True
        else:
            assert key in self.keys(), "filter_by_value():'%s' not found in this expression object" % key
            its_a_condition = False
            
        new_expn = []
        
        conv_dict = {"gt": ">", "lt": "<", "gte": ">=", "lte": " <=", "equal": "=="}
        if evaluator in ("gt", "lt", "gte", "lte", "equal"):
            evaluator = conv_dict[evaluator]
        
        if its_a_condition:
            cond_index = self._conditions.index(key)
            for item in self.linearData:
                if eval("%s %s %s" % (item["conditions"][cond_index], evaluator, value)):
                    new_expn.append(item)
        else: # filter on a normal key.
            for item in self.linearData:
                if eval("%s %s %s" % (item[key], evaluator, value)):
                    new_expn.append(item)
        
        ret = expression(loadable_list=new_expn, cond_names=self._conditions, silent=True)
        
        config.log.info("filter_by_value(): Filtered expression for ['%s' %s %s], found: %s" % (key, evaluator, value, len(ret)))
        return(ret)     

    def filter_low_expressed(self, min_expression, number_of_conditions):
        """
        **Purpose**
            filter genes by a minimum_expression value in at least number_of_conditions
            
            Basically the R command:
            
            keep <- rowSums(expression_table > min_expression) >= number_of_conditions
            
        **Arguments**
            min_expression (Required)
                The minimum expression value required to pass the test
                
            number_of_conditions (Required)
                The number of conditions that must be greater than min_expression
        
        **Results**
            Returns a new genelist 
        """
        
        newl = self.shallowcopy()
        newl.linearData = []
        
        for item in self.linearData:
            if sum([int(i > min_expression) for i in item["conditions"]]) >= number_of_conditions: # passed
                newl.linearData.append(item.copy())
        
        newl._optimiseData()
        config.log.info("filter_low_expression(): removed %s items, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def drawDotplot(self, x_condition_name, y_condition_name, filename=None, genelist=None, key=None, 
        label=False, **kargs):
        """
        Deprecated alias for scatter()
        """
        return(self.scatter(x_condition_name, y_condition_name, filename=None, genelist=None, key=None, 
            label=False, **kargs))

    def scatter(self, x_condition_name, y_condition_name, filename=None, genelist=None, key=None, 
        label=False, **kargs):
        """
        **Purpose**
            draw an X/Y dot plot or scatter plot, get R^2 correlation etc.

        **Arguments**
            x_condition_name = the name of the er... X condition
            y_condition_name = the name of the er... Y condition
            
            genelist (Optional)
                If you send a genelist and a key then these items will be emphasised on the dotplot.
                
            key (Optional, Required if 'genelist' used)
                The key to match between the expression data and the genelist.
                
            label (Optional, default=False)
                If genelist and key are set, then you can label these spots with the label from
                genelist[key].
                If you want to label all of the items in the list, then just send the original genelist
                and a key to use and all of the spots will be labelled.

            do_best_fit_line (Optional, default=False)
                draw a best fit line on the scatter

            print_correlation (Optional, default=None)
                You have to spectify the type of correlation to print on the graph.
                valid are:
                    r = R (Correlation coefficient)
                    r2 = R^2.
                You need to also set do_best_fit_line=True for this to work
                
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
        
        if not "xlabel" in kargs:
            kargs["xlabel"] = x_condition_name   
        if not "ylabel" in kargs:
            kargs["ylabel"] = y_condition_name
        
        if genelist and key:
            matches = genelist.map(genelist=self, key=key) # make sure resulting object is array
            tx = matches.getDataForCondition(x_condition_name)
            ty = matches.getDataForCondition(y_condition_name)
            if label:
                kargs["spot_labels"] = matches[key]
            
            real_filename = self.draw.nice_scatter(x_data, y_data, filename, spots=(tx, ty), **kargs)
        else:
            real_filename = self.draw.nice_scatter(x_data, y_data, filename, **kargs)
        

        config.log.info("scatter(): Saved '%s'" % real_filename)
        return(True)

    def boxplot(self, filename=None, **kargs):
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
        assert filename, "must provide a filename"

        data = self.serialisedArrayDataList

        if "log" in kargs and kargs["log"]:
            data = self.__log_transform_data(data, log=kargs["log"])

        # do plot
        actual_filename = self.draw.boxplot(data=data, filename=filename, labels=self.getConditionNames(), **kargs)

        config.log.info("boxplot(): Saved %s" % actual_filename)
        return(actual_filename)

    def multi_line(self, filename=None, alpha=None, **kargs):
        """
        **Purpose**
            Draw each gene as an alpha blended line. Ideal for visualising continuous data, for example time courses.
        
        **Arguments**
            filename (Required)
                The filename of the image to save.
                
            alpha (Optional, default=guess a value)
                The alpha blending value. By default glbase will attempt to guess a suitable 
                alpha value. But this may be too light/dark for your tastes. You can set the value here, from
                0.0 (completely transparent) to 1.0 (completely opaque)
                
        **Return** 
            None 
        """
        assert filename, "You must specify a filename"
        
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        # Guess an alpha value:
        if not alpha:
            alpha = 8.0 / len(self) 
        
        for item in self:
            ax.plot(item["conditions"], color="grey", alpha=alpha)
            
        self.draw.do_common_args(ax, **kargs)
        realfilename = self.draw.savefigure(fig, filename)
        config.log.info("multi_line(): Saved '%s' with alpha=%.5f" % (realfilename, alpha))
        return(None)

    def drawBoxplot(self, filename=None, **kargs):
        """
        Deprecated alias for boxplot()
        """
        return(self.boxplot(filename=filename, **kargs))

    def log(self, base=math.e, pad=0.00001):
        """
        **Purpose**
            log transform the data

            NOTE: THis is one of the few IN-PLACE glbase commands.

        **Arguments**
            base (Optional, default=math.e)
                the base for the log transform.
            
            pad (Optional, default=1e-6)
                value to pad all values by to log(0) errors.
        
        **Returns**
            None
        """
        do_log = False

        if base == math.e:
            do_log = math.e
        elif isinstance(base, bool):
            do_log = math.e
        elif isinstance(base, int):
            do_log = base
        else:
            do_log = False
            
        for item in self:
            item["conditions"] = [math.log(v+pad, do_log) for v in item["conditions"]]
        self._optimiseData()
        return(None)       

    def unlog(self, base=None, adjuster=0.00001):
        """
        **Purpose**
            retrun the raw data to the unlogged form.  YOU MUST PROVIDE THE CORRECT BASE

            NOTE: THis is one of the few IN-PLACE glbase commands.

        **Arguments**
            base (Required)
                the base for the log transform.
        
        **Returns**
            None
        """          
        for item in self:
            item["conditions"] = [base**(v+adjuster) for v in item["conditions"]]
        self._optimiseData()
        return(None)       

    def mult(self, number=None):
        """
        **Purpose**
            multiply the data by some number
            
            This method came about as I had a reason to multiply RNA-seq data to get it above zero
            Particularly concerning the TPM measure.

            NOTE: This is one of the few IN-PLACE glbase commands.

        **Arguments**
            number (Required)
                the number to multiply by
        
        **Returns**
            None
        """
        assert number, "must provide a number"
           
        for item in self.linearData:
            item["conditions"] = [v*number for v in item["conditions"]]
        self._optimiseData()
        return(None)   

    def add(self, number=None):
        """
        **Purpose**
            add number to the expression data
            
            NOTE: This is one of the few IN-PLACE glbase commands.

        **Arguments**
            number (Required)
                the number to add to the expression data
        
        **Returns**
            None
        """
        assert number, "must provide a number"
           
        for item in self.linearData:
            item["conditions"] = [v+number for v in item["conditions"]]
        self._optimiseData()
        return(None)   

    def abssub(self, number=None):
        """
        **Purpose**
            subtract number from all of the expression data, moving the value towards 0 with a limit of 0.
            
            This will perform the following logic:
            
            if expression_value > 0:
                expression_value -= 1
                if expression_value < 0:
                    expression_value = 0
                    
            elif expression_value < 0:
                expression_value += 1
                if expression_value > 0:
                    expression_value = 0
                
            
            NOTE: This is one of the few IN-PLACE glbase commands.

        **Arguments**
            number (Required)
                the number to add to the expression data
        
        **Returns**
            None
        """
        assert number, "must provide a number"
           
        for item in self.linearData:
            newcond = []
            for expression_value in item["conditions"]:
                if expression_value >= 0:
                    expression_value -= 1
                    if expression_value < 0:
                        expression_value = 0
                    
                elif expression_value <= 0:
                    expression_value += 1
                    if expression_value > 0:
                        expression_value = 0
                newcond.append(expression_value)
            item["conditions"] = newcond
        self._optimiseData()
        return(None)   

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
        config.log.info("curves(): Saved '%s'" % filename)
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

        newl = self.deepcopy()
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
        self._conditions.append(condition_name)
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
            If the expression is below 'minimum_value' in at least
            'number_to_pass' conditions

        **Arguments**

        **Returns**
            A new expression object with probes that fail to pass removed.
        """
        # kargs tidier to go here.

        newl = self.shallowcopy()
        newl.linearData = []

        for probe in self.linearData:
            s = sum([i > minimum_value for i in probe["conditions"]])
            if s > number_to_pass:
                newl.append(probe)
        newl._optimiseData()
        config.log.info("remove_by_expression(): %s items < %s in %s conditions" % (len(self) - len(newl), minimum_value, number_to_pass))
        return(newl)

    def cumulative_distributions(self, genelists=None, filename=None, key=None, **kargs):
        """
        **Purpose**

            draw a set of cumulative distributions, based on a selection
            of genelist-like objects that can be mapped to
            this expression object using 'key'

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
            
            You probably want barh_single_item() in preference over this method
            
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
                The input comes in two forms::
                    
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
        
            Other arguments understood by the figure
                title, xlabel, ylabel
                
                and 'cols' - where you can send a list of colours for the bars. The list
                must be as long as the data, or it will throw an error.
        
        **Example**
        
            Okay, this is pretty complicated (It's because the expression system has only a rudimentary 
            understanding of error values, I may change this in the future)::
            
                x = expression(...)
                
                # I don't have any error values
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
            
        empty_array = [0.0001 for x in self._conditions]
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
                labels=labels, cond_names=self._conditions,
                **kargs) 
        
    def hist(self, filename=None, range=None, suppress_zeros=True, log=None, kde=False, 
        covariance=0.2, bins=50, **kargs):
        """
        **Purpose**
            Draw a normal histogram of the expression values
        
        **Arguments**
            filename (Required)
                the filename to save the resulting image to.
                
            covariance (Optional, default=0.2)
                undocumented wierdness (Sorry)
                
            range (Optional, default=(min(expression), max(expression)))
                This is the maximum value to build the KDE over. (ie. 0 ... mmax).
                You probably really want to chage this value, as small outliers will distort the
                distribution massively.
                Play around with the value until it gets you a nice normal-like distribution.
                
            suppress_zeros (Optional, default=False)
                microarray expression-data don't generally include zero values for expression. 
                However, RNA-seq does include zero values for some transcripts.
                In some instances you may want to suppress these values.
                
                If this is set to True expression values of zero are not included in the 
                histogram.
                
            log (Optional, default=False)
                log transform the data.
                At the moment only log2 is supported.
        
            kde (Optional, default=False)
                use kernel density estimation to smooth the data
                
            covariance (Optional, default=0.2)
                Value for KDE smoothing.
                
            bins (Optional, default=50)
                number of bins to use. (Also affects the resolution of the kde smoothing)
            
        **Returns**

        """
        assert filename, "Need filename, leh!"
        
        binned_data = {}
        
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        all_data = {}
            
        for k in self._conditions:
            expn_values = self.getDataForCondition(k)
            if not range: # sample a range if none specified
                range = (min(expn_values), max(expn_values))
            
            if suppress_zeros:
                expn_values = [v for v in expn_values if int(v*1000000000) != 0]
            
            expn_values = numpy.array(expn_values)
            
            if log:
                expn_values = numpy.log2(expn_values)
            
            if kde:
                expn_values = utils.kde(expn_values, range=range, covariance=covariance, bins=bins)
                ax.plot(expn_values, label=k)
            else:
                ax.hist(expn_values, bins=bins, range=range, normed=True, histtype='step', label=k)
            all_data[k] = expn_values
            
        ax.legend(ncol=1)
                    
        self.draw.do_common_args(ax, **kargs)
        real_filename = self.draw.savefigure(fig, filename)
        
        config.log.info("hist(): Saved '%s'" % real_filename)
        return({"data": all_data, "labels": self._conditions})
        
    def print_stats(self, log=None):
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

        print "All Array:"
        print "\tmin", self.min()
        print "\tmax", self.max()
        print "\tmean", self.mean()
        print "\tstddev", self.stdev()
            
        for k in self._conditions:
            expn_values = self.getDataForCondition(k)
            if log:
                expn_values = [v for v in expn_values if int(v*10000) > 0]
                expn_values = numpy.log2(expn_values)
            
            aaverage = numpy.average(expn_values)
            amedian = numpy.median(expn_values)
            astd = numpy.std(expn_values)
            amin = numpy.min(expn_values)
            amax = numpy.max(expn_values)
            
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
            "conditions": self._conditions})

    def stats(self, log=None):
        """
        **Purpose**
            Get the min, max, mean, median and stddev for each condition in the expression data
        
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
        
        for k in self._conditions:
            expn_values = self.getDataForCondition(k)
            if log:
                expn_values = [v for v in expn_values if int(v*10000) > 0]
                expn_values = numpy.log2(expn_values)
            
            aaverage = numpy.average(expn_values)
            amedian = numpy.median(expn_values)
            astd = numpy.std(expn_values)
            amin = self.min()
            amax = self.max()
            
            means.append(aaverage)
            medians.append(amedian)
            stds.append(astd)
            mins.append(amin)
            maxs.append(amax)
            
        return({"mean": means, "median": medians, "stdev": stds, "min": mins, "max": maxs, 
            "conditions": self._conditions})
        
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

    def ttest(self, condition1, condition2, equal_var=False):
        """
        **Purpose**
            perform a t-test between condition1 and condition2.
            
            Assumes the data is normally distributed. The t-test is only valid if the 
            data is normally distributed. 
            
            This performs and independent t-test of the samples. And assumes the
            two samples are independent. An appropriate use of this test would be 
            to compare between a set of genes across two different biological samples.
            
            See Scipy.stats.ttest_ind() for more details
        
        **Arguments**
            condition1
                the name of condition 1 to use
                
            condition2
                the name of condition 2 to use
                
            equal_var (Optional, default=False)
                Generally the variation in samples is not expected to be the same for a typical
                set of expression data. 
                
                So, set this to False and perform a 'Welch's t-test instead'. You can set it to True
                and perform a Student's independent t-test which assumes sample variances are the same.
                
                See:
                http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
                
                For more details.
                
        **Returns**
            The t-statistic and the p-value
        """
        assert condition1 in self._conditions, "condition1 '%s' not on this array" % condition1
        assert condition2 in self._conditions, "condition2 '%s' not on this array" % condition2
        
        data1 = self.getDataForCondition(condition1)
        data2 = self.getDataForCondition(condition2)
        
        return(ttest_ind(data1, data2, equal_var=equal_var))

    def pearsonr(self, condition1, condition2):
        """
        **Purpose**
            perform a pearson correlation test between condition1 and condition2.
            
            See Scipy.stats.pearsonr() for more details
        
        **Arguments**
            condition1
                the name of condition 1 to use
                
            condition2
                the name of condition 2 to use
                
        **Returns**
            The pearson correlation co-efficient and the p-value
        """
        assert condition1 in self._conditions, "condition1 '%s' not on this array" % condition1
        assert condition2 in self._conditions, "condition2 '%s' not on this array" % condition2
        
        data1 = self.getDataForCondition(condition1)
        data2 = self.getDataForCondition(condition2)
        
        t, p = pearsonr(data1, data2)
        return(t, p)
        
    def kruskal(self, list_of_conditions):
        """
        **Purpose**
            perform a Kruskal-Wallis H-test for independent samples for several conditions.
                                 
            See Scipy.stats.kruskal() for more details
        
        **Arguments**
            list_of_conditions (Required)
                a list of condition names in this expression object to test.
                
        **Returns**
            The kruskal correlation co-efficient and the p-value
        """       
        datas = [self.getDataForCondition(c) for c in list_of_conditions]
        
        t, p = kruskal(*datas)
        return(t, p)

    def mannwhitneyu(self, condition1, condition2):
        """
        **Purpose**
            perform a Mann Whitney U-test for independent samples for two conditions.
                                 
            See Scipy.stats.mannwhitneyu() for more details
        
        **Arguments**
            condition1
                the name of condition 1 to use
                
            condition2
                the name of condition 2 to use
                
        **Returns**
            The u-statistic and one-sided p-value
        """       
        return(self.__unified_stats("mannwhitneyu", condition1, condition2))

    def wilcoxon(self, condition1, condition2):
        """
        **Purpose**
            perform a 'Wilcoxon signed-rank' test between condition1 and condition2.
                                 
            See Scipy.stats.wilcoxon() for more details
        
        **Arguments**
            condition1
                the name of condition 1 to use
                
            condition2
                the name of condition 2 to use
                
        **Returns**
            The wilcoxon z-statistic and p-value
        """
        return(self.__unified_stats("wilcoxon", condition1, condition2))

    def __unified_stats(self, test, condition1, condition2):
        assert condition1 in self._conditions, "condition1 '%s' not on this array" % condition1
        assert condition2 in self._conditions, "condition2 '%s' not on this array" % condition2
        
        data1 = self.getDataForCondition(condition1)
        data2 = self.getDataForCondition(condition2)
        
        test_funcs = {"wilcoxon": stats.wilcoxon,
            "mannwhitneyu": stats.mannwhitneyu,
            }
        
        return(test_funcs[test](data1, data2))            
        
    def tree(self, mode="conditions", filename=None, row_names=None, 
        cluster_mode="euclidean", color_threshold=None, label_size=7, **kargs):
        """
        **Purpose**
            Draw a hierarchical clustered tree of either the 'conditions' or 'rows'
            
        **Arguments**
            filename (Required)
                filename to save an image to.
            
            mode (Optional, default="conditions")
                cluster either the "conditions" or the "rows". (rows are usually genes)
                
            row_names (Optional, default=None)
                A key to use for the row_names, if mode == "row"
                
            cluster_mode (Optional, default="euclidean")
                A metric to cluster the data by.
                
            color_threshold (Optional, default=None)
                By default tree() uses the Scipy/MATLAB default colours to signify
                links with similar distances. Set to -1 to change the tree to all blue.
                
            label_size (Optional, default=7)
                The size of the text attached to the leaf labels.
                
        **Returns**
            The Tree in a dictionary {"tree": tree_data, "a": <dendrogram_data>}
        """
        valid_modes = ("conditions", "rows")
        
        assert mode in valid_modes, "'%s' not a valid mode" % mode
        
        if not "size" in kargs: # resize if not specified
            kargs["size"] = (3,6)
        
        fig = self.draw.getfigure(**kargs)
        
        data = numpy.array(self.serialisedArrayDataList) 
            
        if mode == "conditions": # Use the condition names for rows:
            row_names = self._conditions
        elif mode == "rows":
            data = data.T
            if row_names:
                row_names = self[row_names]

        ax = fig.add_subplot(111)
        ax.set_position([0.01, 0.01, 0.4, 0.98])

        # from scipy;
        # generate the dendrogram
        Y = pdist(data, metric=cluster_mode)
        Z = linkage(Y, 'complete')
        a = dendrogram(Z, orientation='right', color_threshold=0.0, labels=row_names)
        
        ax.set_frame_on(False)
        ax.set_xticklabels("")
        [item.set_markeredgewidth(0.0) for item in ax.xaxis.get_ticklines()]
        #[item.set_markeredgewidth(0.0) for item in ax.yaxis.get_ticklines()]

        # Use the tree to reorder the data.
        row_names = a["ivl"]
        [t.set_fontsize(label_size) for t in ax.get_yticklabels()]
                
        #self.draw.do_common_args(ax, **kargs) # broken!?
        if filename: 
            real_filename = self.draw.savefigure(fig, filename)
            config.log.info("tree(): Saved '%s'" % real_filename)
        return({"tree": Z, "a": a, "Z": Z})
        
    def gene_curve(self, key=None, values=None, filename=None, moving_average=None,
        **kargs):
        """
        **Purpose**
            Draw all of the genes specified on a single axis.
            
        **Arguments**
            key (Required)
                The key to use to match to items found in values
                
            values (Required)
                The values to match, should be a list or other iterable
                
            filename (Required)
                the filename to save to
                
            moving_average (Optional, default=False)
                use a moving average to smooth the gene curves.
                If set then should be an integer specifying the number of neighboring bins to use
        
        **Returns**
            A file in filename and None
        """
        assert key, "Must specify key"
        assert values, "must specify values"
        assert filename, "must specify filename to save to"
        
        if not "aspect" in kargs:
            kargs["aspect"] = "wide"
        
        keeps = frozenset(values)
        
        x_data = numpy.arange(len(self._conditions))
        
        data = {}
        for i in self.linearData:
            if i[key] in keeps:
                data[i[key]] = i["conditions"]
        
        fig = self.draw.getfigure(**kargs)
        
        ax = fig.add_subplot(111)
        
        for k in data:
            if moving_average:
                x_data, y_data = utils.movingAverage(data[k], window=moving_average)
                ax.plot(x_data, y_data, label=k)
            else:
                ax.plot(x_data, data[k], label=k)
        
        ax.legend()
        self.draw.do_common_args(ax, **kargs)
        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("gene_curve(): Saved '%s'" % actual_filename)
        
    def correlation_heatmap(self, filename=None, mode="r2", aspect="square", bracket=(0,1), **kargs):
        """
        **Purpose**
            Plot a heatmap of the (R, R^2, Pearson or Spearman) correlation for all pairs of
            samples in the expression object
            
        **Arguments**
            filename (Required)
                the filename to save the image to.
                
            mode (Optional, default="r2")
                by default the R (Coefficient of determination) is squared. Set to 'r' for 
                Coefficient of determination value.
                use 'pearson' for a Pearson correlation score and 'spearman' for a 
                Spearman-ranked score.   
                
            bracket (optional, default=(0,1))
                bracket the heatmap by these values.
                
            Other heatmap options that will work:
                heat_hei
                heat_wid
                ...
                
        **Returns**
            A dict containing:
                "data": 2D numpy array containing the correlation scores (depending upon the mode)
                "p": 2D numpy array of p-values.
                "labels": the labels along the top and bottom of the array, sorted according to the 
                    clustering (if any)
                
            Note that the arrays are not clustered.
        """
        assert filename, "You must specify a filename"
        
        arr = numpy.zeros((len(self._conditions), len(self._conditions)))
        ps = numpy.zeros((len(self._conditions), len(self._conditions)))
        
        p = progressbar(len(self._conditions))
        for ic1, c1 in enumerate(self._conditions):
            for ic2, c2 in enumerate(self._conditions):
                if c1 == c2:
                    arr[ic1, ic2] = 1.0
                else:
                    x = self.getDataForCondition(c1)
                    y = self.getDataForCondition(c2)
                    if mode in ("r", "r2"):
                        arbr = scipy.polyfit(x, y, 1)
                        xr = scipy.polyval(arbr, x)
                        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
                        ps[ic1, ic2] = p_value
                        if mode == "r":
                            arr[ic1, ic2] = r_value
                        elif mode == "r2":
                            arr[ic1, ic2] = r_value**2
                            
                    elif mode == "pearson":
                        arr[ic1,ic2], ps[ic1, ic2] = scipy.stats.pearsonr(x, y)
                    elif mode == "spearman":
                        arr[ic1,ic2], ps[ic1, ic2] = scipy.stats.spearmanr(x, y)
            p.update(ic1)
        
        square = True
        aspect = "square"
        if "heat_hei" in kargs or "heat_wid" in kargs:
            square=False
                
        results = self.draw.heatmap(filename=filename, data=arr, square=square,
            bracket=bracket, aspect="square", row_names=self._conditions, col_names=self._conditions, 
            colbar_label="correlation (%s)" % mode, **kargs)
        config.log.info("correlation_heatmap(): Saved '%s'" % results["real_filename"])
        return({"data": results["reordered_data"], "labels": results["reordered_cols"]})

    def closest_correlate(self, target_condition, number=5, cut_off=0.7, method="r2", **kargs):
        """
        
        **Purpose**
            Find the closest correlating samples to <target_condition>
            
        **Arguments**
            target_condition (Required)
                The name of the target condition to find the closest correlates for
                
            number (Optional, default=5)
                report <number> closest correlates
            
            cut_off (Optional, default=0.7)
                Do not report a correlate if it falls below the <cut_off> value
                
            method (Optional, default="r2")
                The method to use for the correlation, must be one of:
                "r", "r2", "spearman", "pearson"
                
        **Returns**
            A list of dictionaries containing the closest correlations, their rank
            and their correlation score
        
        """
        assert target_condition in self._conditions, "'%s' condition was not found" % target_condition
        assert method in ("r", "r2", "pearson", "spearman"), "method '%s' not recognised" % method
        
        # get all the correlations:
        
        res = []
        x = self.getDataForCondition(target_condition)
        
        for ic1, c1 in enumerate(self._conditions):
            if c1 != target_condition:
                y = self.getDataForCondition(c1)
                if method in ("r", "r2"):
                    arbr = scipy.polyfit(x, y, 1)
                    xr = scipy.polyval(arbr, x)
                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
                    if method == "r":
                        correlation = r_value
                    elif method == "r2":
                        correlation = r_value**2
                    
                elif method == "pearson":
                    correlation, p_value = scipy.stats.pearsonr(x, y)
                elif method == "spearman":
                    correlation, p_value = scipy.stats.spearmanr(x, y)
            
                if correlation > cut_off:
                    res.append({"name": c1, "correlation": correlation})
                
        # sort the dict by correlation:
        res = sorted(res, key=itemgetter("correlation"))
        res.reverse()
        return(res)

    def cut(self, function):
        """
        **Purpose**
            Perform an action (any function(), often a lambda) on each item 
            and return the resulting expression object.
            
            This function does (essentially) this::
            
            newlist = []
            for item in self:
                if function(item):
                    newlist.append(item)
            return(newlist) 
                
        **Arguments**
            function (Required)
                A Python method or lambda, that takes one argument, which is each row
                of the item. Must return True or False
                
        **Returns**
            A new expression object.
        """
        
        assert function, "you must supply a falid function"
        assert function({"name": None}), "function appears to be invalid"
        
        newl = []
        for item in self.linearData:
            print item
            print function(item)
            if function(item):
                newl.append(item)
        
        if not newl:
            config.log.info("cut(): result of cut() is empty!")
            return(None)
        
        cc = self.shallowcopy() # going to replace linearData.
        cc.linearData = newl
        cc._optimiseData()
        
        return(cc)
        
    def barh_single_item(self, key=None, value=None, filename=None, tree=None, 
        plot_mean=True, plot_stdev=True, fold_change=False, **kargs):
        """
        **Purpose**
            Plot a horizontal bar for a single item (gene?), for all conditions.
            
            Uses an 'err' key for error bars if present
            
            Note that this will find only the first item in the list that matches.
            
        **Arguments**
            key (Required)
                The key to search for the item
                
            value (Required)
                the value to search for
                
            filename (Required)
                the filename to save the plot to
                
            tree (Optional, default=None)
                A tree used to order the conditions by. Should be the output from tree()
                and should only be clustering the conditions
            
            plot_mean (Optional, default=True)
                plot a grey line showing the mean of all of the expression values
                
            plot_stdev (Optional, default=True)
                plot a blue line at 1x stdev and a red line at 2x stdev
            
            fold_change (Optional, default=False)
                by default barh_single_itme expects absolute levels of expression, but if you want to
                provide fold-change data (which centres around 0 then setting this to True will
                modify the plot and make it more friendly for fold-change plost
                (re-bracket, a line at 0, and 2 and 4 fold-change (assuming data is
                log2(fold-change)), plot_mean and plot_stdev=False).
        
            Typical arguments for plots are supported:
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits (Note that barh_single_item will clamp on [0, max(data)+10%] by default.
                ylims - y-axis limits
                xticklabel_fontsize - x tick labels fontsizes
                yticklabel_fontsize - y tick labels fontsizes
                
        
        **Returns** 
            The item used to draw the barchart    
        """
        assert key, "you must specify a key"
        assert value, "you must specify a value"
        assert filename, "you must specify a filename"
        
        item = self._findDataByKeyLazy(key, value)
        
        if not item:
            config.log.warning("'%s:%s' not found in this list, not saving" % (key, value))
            return(None)
        
        err = None  
        # re-order "conditions" by tree, if present:
        if tree:
            data = []
            conds = []
            if "err" in item:
                err = []
            order = tree["a"]["ivl"]

            # resort the data by order;
            for label in order:
                indx = self._conditions.index(label)
                data.append(item["conditions"][indx])
                conds.append(label)
                if "err" in item:
                    err.append(item["err"][indx])

        else:
            data = list(item["conditions"])
            conds = list(self._conditions)
            if "err" in item:
                err = list(item["err"])
         
        fig = self.draw.getfigure(aspect="long", **kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.46, 0.04, 0.5, 0.92])
        
        y = numpy.arange(len(data))
        if err:
            # error bars'd graphs look better with black on grey
            ax.barh(y, data, xerr=err, ecolor="black", ec="none", fc="grey", height=0.5)
        else:
            # no error bars, solid black is better
            ax.barh(y, data, ec="none", fc="black", height=0.5)
        ax.set_yticklabels(conds)
        ax.set_yticks(y+0.25) # 0.25 should be half of the height of the bars so that text aligns with the bar
        ax.set_ylim([0, len(data)])
        [item.set_markeredgewidth(0.0) for item in ax.yaxis.get_ticklines()]
        
        if fold_change:
            ax.axvline(-1, color="red", ls=":")
            ax.axvline(1, color="red", ls=":")
            ax.axvline(-2, color="orange", ls=":")
            ax.axvline(2, color="orange", ls=":")
            ax.axvline(-3, color="yellow", ls=":")
            ax.axvline(3, color="yellow", ls=":")
        else:
        
            if plot_mean:
                ax.axvline(numpy.average(data), color="grey", ls=":")
                ax.text(numpy.average(data)+0.1, len(data)+0.02, "m", size=6)
        
            if plot_stdev:
                mean = numpy.average(data)
                std = numpy.std(data)
                ax.axvline(mean + std, color="blue", ls=":")
                ax.axvline(mean - std, color="blue", ls=":")
                ax.axvline(mean + (std*2.0), color="red", ls=":")
                ax.text(mean + std+0.1, len(data)+0.02, "1x", size=6)
                ax.text(mean + (std*2.0)+0.1, len(data)+0.02, "2x", size=6)

        if not "xlims" in kargs:
            if fold_change:
                mx = max(data)+(max(data)*0.1)
                if mx < 2.0:
                    mx = 2.0
                ax.set_xlim([-mx, mx])
            else:
                ax.set_xlim([0, max(data)+(max(data)*0.1)])
        
        if fold_change:
            ax.axvline(0, color="grey", ls="-")
            
        if not "yticklabel_fontsize" in kargs:
            kargs["yticklabel_fontsize"] = 10
            
        if not "title" in kargs:
            kargs["title"] = value
               
        self.draw.do_common_args(ax, **kargs)
        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("barh_single_item(): Saved '%s'" % actual_filename)
        return(item)
        
    def get_pca(self, rowwise=False, label_key=None, **kargs):
        """
        **Purpose**
            Get the pca object for this expression data. For performing principal 
            component analysis (PCA) on the expression data set.
            
            A typical workflow for PCA of expression data would go something like:
            
            expn = expression(filename="...", ...)
            pca = expn.get_pca()
            print pca.max() # maximum number of PCs.
            pca.loading(filename="...") # PC loading bar chart.
            pca.plot(1,2, filename="...") # PC 
            
        **Arguments**
            rowwise (Optional, default=False)
                Perform PCA on the rows (probably genes), instead of on the conditions (the default)
                
            label_key (Optional, Required if rowwise=True)
                The key to use in the genelist to use as a label on plots, etc. 
                
                Ignored if rowwise=False
            
        **Returns**
            A 'pca' object.
            See the documentation for pca for more details. 
            
            You can get the pca object again in expn.pca
        """
        self.pca = pca(self, rowwise=rowwise, label_key=label_key)
        return(self.pca)
        
    def fc_scatter(self, cond1, cond2, filename=None, plot_diagonals=True, zoom_bracket=[-3,3], 
        label_key="name", **kargs):
        """
        **Purpose**
            Draw a scatter plot of the fold-changes between two conditions
            
            There is not so much difference between this and scatter()
            The major advantage is it is set up to label particular genes, genes moving away
            from the diagonals and has more sensible defaults for a scatter plot based on fold-change.
            
            You could certainly emulate this using the normal scatter(). But think of this as a nice helper 
            function.
            
            NOTE: This does NOT fold-change transform your data and assumes
            you have already done this. It also does not check if the data is fold-change data.
            
            use something like norm_multi_fc() or normaliseToCondition() before using this function.
            
            By default it will plot two plots, one zoomed in and one zoomed out to show all of the data
            
        **Arguments**
            cond1 (Required)
                the name of condition1
                
            cond2 (Required)
                the name of condition2
                
            filename (Required)
                the filename to save the scatter plot to
                
            plot_diagonals (Optional, default=True)
                plot guidelines showing the diagonals
                
            zoom_bracket (Optional, default=-[3,3])
                the bracket to use for the zoomed in plot
                
            label_key (Optional, default="name")
        
        **Returns**
            None
        """
        assert filename, "fc_scatter(): you must specify a filename"
        
        # get the data
        c1d = self.getDataForCondition(cond1)
        c2d = self.getDataForCondition(cond2)
        names = self[label_key]
        
        # load the data for plotting
        pt_x = []
        pt_y = []
        labs = []
        cols = []
        for index, x in enumerate(c1d):
            pt_x.append(c1d[index])
            pt_y.append(c2d[index])
            cols.append("blue")
                
            # do the colouring here if away from the diagonal
            if label_key:
                labs.append(names[index])
            """
                    # get the distance from the diagonal:
            dx = (ex[l1] - ex[l2]) # i.e. x - y
            dy = (ex[l2] - ex[l1]) # i.e. y - x

            d = dx * dx # actually the squared distance
            print dx, dy, d    
            if d > 0.05 or d < -0.05:
                cols.append("red")
                labs.append(g["name"])
                print "%s\t%.2f" % (g["name"], d)
                both.append({"ensg": g["ensg"]})
            else:
                cols.append("blue")
                labs.append("")
            """
        
        # plot
        if not "size" in kargs:
            kargs["size"] = (13,6)
        fig = self.draw.getfigure(**kargs)
        
        ax = fig.add_subplot(121)    
        ax.scatter(pt_x, pt_y, alpha=0.2, color=cols)
        for i in xrange(len(labs)):
            ax.text(pt_x[i], pt_y[i], labs[i], size=5, ha="center")

        # Diagonal slopes:
        ax.plot([5, -5], [5, -5], ":", color="grey")
        ax.plot([-5, 5], [5,-5], ":", color="grey")

        ax.axvline(0, color="grey", ls="-.")
        ax.axhline(0, color="grey", ls="-.")
        ax.set_xlim(zoom_bracket)
        ax.set_ylim(zoom_bracket)
        ax.set_xlabel(cond1)
        ax.set_ylabel(cond2)

        ax = fig.add_subplot(122)    
        ax.scatter(pt_x, pt_y, alpha=0.2, color=cols)
        #for i in xrange(len(labs)):
        #    ax.text(pt_x[i], pt_y[i], labs[i], size=6, ha="center")

        # Diagonal slopes:
        if plot_diagonals:
            ax.plot([5, -5], [5, -5], ":", color="grey")
            ax.plot([-5, 5], [5,-5], ":", color="grey")

        minmax = max([abs(min(c1d)), max(c1d), abs(min(c2d)), max(c2d)]) 

        ax.axvline(0, color="grey", ls="-.")
        ax.axhline(0, color="grey", ls="-.")
        ax.set_xlim([-minmax,minmax])
        ax.set_ylim([-minmax,minmax])
        ax.set_xlabel(cond1)
        ax.set_ylabel(cond2)

        # best fit:
        """
        (ar, br) = polyfit(pt_x, pt_y, 1)
        print ar, br
        #xr = polyval([ar,br], pt_x)
        slope, intercept, r_value, p_value, std_err = linregress(pt_x,pt_y)

        print "m, x, r2:", slope, intercept, r_value*r_value

        mx = [min(pt_x), max(pt_x)]
        my = [(slope * min(pt_x)) + intercept, (slope * max(pt_x)) + intercept]

        ax.plot(mx, my, "r.-")
        """

        actual_filename = self.draw.savefigure(fig, filename)     
        config.log.info("fc_scatter(): Saved '%s'" % actual_filename)
        return(None)
        
    def kmeans(self, filename=None, key=None, seeds=None, dowhiten=True, plotseeds=True, **kargs):
        """
        **Purpose**
            perform k-means clustering on the expression data. 
            
            There are two modes to k-means: Either randomly gather 
            patterns from the data or the user can provide the initial patterns to 
            cluster around.
            
            Output a stack of plots and return the groups
            
        **Arguments**
            filename (Optional)
                filename to save the profiles to
            
            key (Optional)
                A key to use to extract the initial values from
                
            seeds (Optional)
                If you send a number:
                    kmeans() will select 0..seeds randomly from the data
                If seeds is a genelist
                    For each item in the genelist, k-means cluster around that 
                    value.
                    
            dowhiten (Optional, default=True)
                'whiten' the data, i.e. standardise the variation. 
                This can be a bad idea if your data is already heavily normalised (i.e.
                variation has already been controlled).
                
            plotseeds (Optional, default=False)
                don't plot the seeds
                    
        **Returns**
            A dictionary containing keys either from the names of the seeds or the random
            group from the seed
        """
        
        # expression table needs to be transformed.
        data = numpy.array(self.serialisedArrayDataList).T
        
        print "data", data.shape
        
        if dowhiten:
            config.log.info("kmeans(): whiten...")
            wt = whiten(data) # features are columns
        else:
            wt = data

        if isinstance(seeds, int):
            seeds = seeds
            seed_names = None
        else: # assume genelist
            assert key, "You must specify a key"
            cents = []
            seed_names = []
            for i in seeds:
                col = self.index(key, i[key])# I need to know the column numbers for each item in seeds: 
                cents.append(self[col]["conditions"])    
                seed_names.append(i[key])
            seeds = numpy.array(cents)
            
            print "seeds", seeds.shape
        
        config.log.info("kmeans(): kmeans...")
        centroids, variance = kmeans(wt, seeds)
        
        print centroids.shape
        print centroids
        
        config.log.info("kmeans(): vq...")
        code, distance = vq(wt, centroids)

        print code

        clusters = {}
        for feature, cluster in enumerate(code):
            if not cluster in clusters:
                clusters[cluster] = []
            clusters[cluster].append(data[feature])
    
        print clusters.keys()
        
        # Guess a suitable arrangement
        sq = math.ceil(math.sqrt(len(clusters)))
        if not "size" in kargs:
            kargs["size"] = (sq*2, sq*2)
        
        fig = self.draw.getfigure(**kargs)

        fig.subplots_adjust(0.01, 0.01, 0.98, 0.98, wspace=0.15, hspace=0.15)

        for i, k in enumerate(clusters):
            ax = fig.add_subplot(sq, sq, i+1)
            ta = numpy.array(clusters[k])
            ax.plot(ta.T, alpha=0.1, color="grey")
                        
            if seed_names:
                # plot he original centroid:
                if plotseeds:
                    ax.plot(cents[i], color="black")
                ax.set_title("%s (%s)" % (seed_names[i], len(clusters[k])), size=6)
            else:
                ax.set_title("Cluster %s (%s)" % (k+1, len(clusters[k])), size=6)
                
            ax.set_xlim([-1,len(self._conditions)])
            #ax.set_xticklabels("")
            #ax.set_yticklabels("")
            
        actual_filename = self.draw.savefigure(fig, filename)     
        config.log.info("kmeans(): Saved '%s'" % actual_filename)