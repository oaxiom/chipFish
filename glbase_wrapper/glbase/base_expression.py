"""

Expression class 

Basic handling for microarray and rna-seq and realtime PCR like data

"""

import sys, os, csv, string, math, copy

from numpy import array, arange, meshgrid, zeros, linspace, mean, object_, std
from array import array as qarray # this should be deprecated later.

import config
from flags import *
from draw import draw
from genelist import genelist
from progress import progressbar
from errors import AssertionError, ArgumentError

class base_expression(genelist):
    def __init__(self, loadable_list=None, filename=None, format=None, expn=None, **kargs):
        """
        **Purpose**
            The base container for expression data.

            Please not though that:
            Expression analysis in glbase is very underpowered and requires normalised
            and easily manipulatable data. Examples for normalisation are
            genespring output, or output from R and PMA, Cufflinks etc...

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
            
            cv_err (Optional)
                Some sort of descriptor for where to get confidence intervals from.
                This should be a tuple, as confidence intervals are 'lo' and 'hi'.

            conditions_names (Optional)
                A list of the condition names (in order) if glbase is not working them 
                out itself.
            
            name (Optional)
                By default expression will use the filename (removing any .txt, .tsv, etc) from
                the ends of the names      
        """
        assert expn, "'expn' argument cannot be empty"
        if not loadable_list:
            assert filename, "no filename to load"
            assert format, "required argument 'format' is missing"
            assert os.path.exists(os.path.realpath(filename)), "'%s' not found" % filename
        else:
            # probably should put some more sanity checking in here.
            assert loadable_list[0], "the list to load does not appear to be a proper list"
        
        if "cv_err" in kargs or "err_up" in kargs or "err_dn" in kargs:
            raise NotImplementedError, "Whoops! I haven't finished expression class - cv_err, err_up and err_dn are not implemented"
            
        valig_args = ["condition_names", "name", "force_tsv"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.__init__, k)

        genelist.__init__(self)
        
        self.filename = filename
        
        self.name = "None"
        if filename:
            self.name = "".join(self.filename.split(".")[:-1])
        if "name" in kargs and kargs["name"]:
            self.name = kargs[k]
        
        if loadable_list:
            self.load_list(loadable_list, expn)
        else:
            # This is a placeholder at the moment,
            # I reload the expn and err values back into the format
            # When you redo this, remember to also redo load_list()
            newf = format
            newf["conditions"] = {"code": expn}
            if "err" in kargs and kargs["err"]:
                newf["err"] = {"code": kargs["err"]}
            elif "cv_err" in kargs and kargs["cv_err"]:
                newf["cv_err"] = kargs["cv_err"]
            
            if "force_tsv" in kargs and kargs["force_tsv"]:
                newf["force_tsv"] = True
            
            self.conditions = [] # Provide a dummy conditions temporarily
            self.loadCSV(filename=filename, format=newf) # no need for error checking here - it's in genelist now.
    
            if "condition_names" in kargs and kargs["condition_names"]:
                self.conditions = kargs["condition_names"]
            else:
                # re-open the file and try to guess the conditions
                # reopen the file to get the condition headers.
                oh = open(filename, "rU")
                if "force_tsv" in format and format["force_tsv"]:
                    reader = csv.reader(oh, dialect=csv.excel_tab)
                if "dialect" in format:
                    reader = csv.reader(oh, dialect=format["dialect"])
                else:
                    reader = csv.reader(oh)
        
                self.conditions = []
                for column in reader:
                    exec "names = %s" % format["conditions"]["code"] # yay, more nice happy arbitrary code execution.
        
                    if names:
                        self.conditions = [str(k) for k in names]
                    break
                oh.close()
                config.log.info("I guessed the condition names as '%s'" % ", ".join(self.conditions))
            
        # coerce the conditions werrs etc to floats
        for i in self:
            i["conditions"] = [float(t) for t in i["conditions"]]
            if "err" in i:
                i["err"] = [float(t) for t in i["err"]]
            if "cv_err" in i:
                i["cv_err"] = [float(t) for t in i["cv_err"]]
        
        self._optimiseData()
        config.log.info("Loaded expression data, found %s items" % len(self))
        
    def __repr__(self):
        return("glbase.expression")

    def _optimiseData(self):
        """
        (Override)
        (Internal)
        Add expression optimisations
        """
        genelist._optimiseData(self) # do the parent optimise.

        # generate a serialised version of the array conditions.

        data = {}
        for array_data in self["conditions"]:
            for index, name in enumerate(self.conditions):
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
        self.serialisedArrayDataList = [self.serialisedArrayDataDict[key] for key in self.conditions]

        return(True)
        
    def saveTSV(self, filename=None, tsv=True, **kargs):
        """
        (Override)
        **Purpose**
            Save the microarray data as a tsv file
            This is a little different from the normal genelist.saveTSV()
            as I want to make certain that the condition data is written in a sensible manner at
            the end of the TSV. 
            I also need to deal with grid like structures etc.
            
            As a general warning, use expression.save() in preference to this. 
            This save is not guaranteed to survive reloading into glbase

        **Arguments**

            filename
                The filename (with a valid path) to save the file to.

        **Returns**

            returns None
            A saved csv file in filename.
        """
        valig_args = ["filename", "tsv", "key_order"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.saveCSV, k)

        assert filename, "you must specify a filename"

        oh = open(os.path.realpath(filename), "w")
        if tsv:
            writer = csv.writer(oh, dialect=csv.excel_tab)
        else:
            writer = csv.writer(oh)

        array_data_keys = ("conditions", "err", "cv_err")

        write_keys = []
        if "key_order" in kargs:
            write_keys = kargs["key_order"]
            # now add in any missing keys to the right side of the list:
            for item in self.keys():
                if item not in write_keys and item not in array_data_keys: # But omit the array_data_keys
                    write_keys.append(item)
        else:
            # just select them all:
            write_keys = self.keys()
            
        title_row = []
        for k in write_keys:
            if k in self.linearData[0]: # sample to see if we have this key
                title_row.append(k) # should mimic below
        writer.writerow(write_keys + self.getConditionNames())

        for data in self.linearData:
            line = []
            for key in write_keys:
                if key in data:
                    line.append(data[key])
            writer.writerow(line + data["conditions"])# conditions go last.
        oh.close()
        
        config.log.info("Saved a csv file to '%s'" % filename)
        return(None)
        
    def sort(self, key):
        """
        This is slightly different from the vanilla genelist's sort - you can pass it the name of
        a condition. Take care to make sure the condition name is not also a valid list key.
        The algorithm searches the genelist before searching the array for your particular condition.

        Also take care with this one: It is one of the few in-place list
        modifiers.

        **Arguments**

        key

            must be a valid key in the genelist or the name of an array condition.

        **Result**

        returns True if succesful.

        returns False if not valid.
        """
        assert (key in self.linearData[0]) or key in self.conditions, "'%s' search key not found in list or array data" % key

        if key in self.linearData[0]:
            return(genelist.sort(self, key)) # use the parents sort.
        else:
            if key in self.conditions:
                name_index = self.conditions.index(key)
                self.linearData = sorted(self.linearData, cmp=lambda x, y: cmp(x["conditions"][name_index],y["conditions"][name_index])) # the original sort() was overridden.
                self._optimiseData()
                return(True)
        return(False)
        
    def getColumns(self, return_keys=None, discard_expression_data=False):
        """
        **Purpose**
            Get and return a new expression-like object containing only the 
            keys specified in return_keys. 
            Note by default the array will also contain the expression data, and the return object
            will be an expression-like object. 
            If 'discard_expression_data=False' then the expression data will be deleted and 
            instead this will return a vanilla genelist object. (It can't be an expression-like object
            anymore as there is no longer any expression data attached to the list).
            
        **Arguments**
        return a new geneList only containing the columns specified in return _keys (a list)
        """
        assert isinstance(return_keys, list), "return_keys must have a list"

        if discard_expression_data:
            newl = self.__copy__()
            newl.linearData = []
            for akey in ["condtions", "err", "cv_err"]: # array data keys to copy over
                if akey in self.linearData[0]:
                    if not akey in return_keys:
                        return_keys.append(akey)
        else:
            newl = genelist() # Discard the expression class and retype as vanilla genelist
            # This may not always be up to date...
            # TODO: It would be better to add some method to transfer generic meta data...
            newl.name = self.name
            newl._history = self._history # a list of historyItem's
            newl.metadata = self.metadata
            newl._history.append("Retyped from '%s' to a genelist" % self.__repr__())

        for item in self.linearData:
            newd = {}
            for key in return_keys:
                newd[key] = item[key]
            newl.linearData.append(newd)
        newl._optimiseData()
        
        self._history.append("Column sliced, for the columns: %s" % return_keys)
        config.log.info("Column sliced, for the columns: %s" % ", ".join(return_keys))
        return(newl)

    def load_list(self, list_to_load, expn=None, name=False):
        """
        **Purpose**
            You've generated your own [{ ... }, { ...}] like list
            (A list of dicts) and you want to either reload it into
            a genelist-like object or load it into an empty genelist.
            This is the method to do that officially.

            This method should be used with great care. Some sanity
            checking is done. But not very much.
            
            This load_list is modified for expression-like genelists. 
            (eg. microarray()). Here you can load keys into conditions based on
            their key names.

        **Arguments**
            list_to_load
                must be a list of dicts.
                
			expn (optional)
				A list of key names to construct the expression data from
				If not specified then it assumes your list already has a correctly formatted
				"conditions" key. 

        **Returns**
            None. This is one of the few IN PLACE methods. and returns
            None.
        """
        assert list_to_load[0], "list_to_load does not appear to be a valid list"

        if expn:
            assert isinstance(expn, list), "'expn' must be a list of keys"
            # Bodge in a new "conditions" key:
            newl = []
            for i in list_to_load:
                new = i.copy()
                new["conditions"] = [i[k] for k in expn]
                for k in expn:
                    del new[k]
                newl.append(new)
            self.conditions = expn
        else:
            newl = list_to_load
                
        # Now call parent with new list
        genelist.load_list(self, newl, name)