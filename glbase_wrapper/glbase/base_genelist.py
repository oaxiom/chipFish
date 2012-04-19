
import copy, cPickle

import config
from helpers import *
from errors import AssertionError, UnRecognisedCSVFormatError, UnrecognisedFileFormatError, ArgumentError

class _base_genelist:
    def __init__(self):
        """
        (Internal)
        This is the base derived class for all genelists.
        It contains methods available to all implementations of genelist.
        """
        self.name = None
        self.linearData = None

    def __repr__(self):
        return("<base genelist class>")

    def __in__(self, key):
        """
        (Override)
        
        Confer: 
        if "key" in genelist:
        """
        return(key in self.keys)

    def __nonzero__(self):
        """
        Fixes:
        if genelist_object:
            Now True
            
        and fixes 
        if genelist: # fails if it is an empty genelist.
        """
        return(len(self) > 0)

    def __copy__(self):
        """
        (Override)
        Confer copy to mean a deep copy as opposed to a shallow copy.
        
        This is required as genelists are compound lists.
        """
        try:
            return(cPickle.loads(cPickle.dumps(self, -1))) # This is 2-3x faster and presumably uses less memory
        except PicklingError:
            return(copy.deepcopy(self)) # Use the generic version

    def __shallowcopy__(self):
        """
        (New)

        Some weird behaviour here, I know, this is so I can still get access to
        the shallow copy mechanism even though 90% of the operations are copies.
        """
        return(copy.copy(self))
        
    def __len__(self):
        """
        (Override)
        get the length of the list
        """
        return(len(self.linearData))
        
    def __int__(self):
        """
        (Override)
        get the length of the list
        NOTE: It's possible this is a bug/feature.
        I don't remove it at the moment as I'm not sure if it is used anywhere.
        
        """
        return(len(self.linearData))
        
    def __iter__(self):
        """
        (Override)
        make the geneList behave like a normal iterator (list)
        """
        for n in self.linearData:
            yield n

    def keys(self):
        """
        return a list of all the valid keys for this geneList
        """
        return([key for key in self.linearData[0]])
     
    def history(self):
        """
        get the origins and history and where this list has been, and all of the operators
        performed on it.
        """
        print self._history

    def _guessDataType(self, value):
        """
        (Internal)

        Take a guess at the most reasonable datatype to store value as.
        returns the resulting data type based on a list of logical cooercions
        (explain as I fail each cooercion).
        Used internally in _loadCSV()
        I expect this will get larger and larger with new datatypes, so it's here as
        as a separate function.

        Datatype coercion preference:
        float > int > location > string
        """
        try: # see if the element is a float()
            if "." in value: # if no decimal point, prefer to save as a int.
                return(float(value))
            else:
                raise ValueError
        except ValueError:
            try: # see if it's actually an int?
                return(int(value))
            except ValueError:
                try: # see if I can cooerce it into a location:
                    return(location(loc=value))
                except (TypeError, IndexError, AttributeError, AssertionError, ValueError): # this is not working, just store it as a string
                    return(str(value))
        return("") # return an empty datatype. 
        # I think it is possible to get here. If the exception at int() or float() returns something other than a 
        # ValueError (Unlikely, Impossible?)

    def _processKey(self, format, column):
        """
        (Internal)
        the inner part of _loadCSV() to determine what to do with the key.
        Better in here too for security.
        """

        d = {}
        for key in format:
            if not (key in ignorekeys): # ignore these tags
                #if not key in d:
                #    d[key] = {}
                    
                if isinstance(format[key], dict) and "code" in format[key]:
                    # a code block insertion goes here - any valid lib and one line python code fragment
                    # store it as a dict with the key "code"
                    d[key] = eval(format[key]["code"]) 
                elif isinstance(format[key], str) and "location" in format[key]:
                    # locations are very common, add support for them out of the box:
                    d[key] = eval(format[key])
                else:
                    d[key] = self._guessDataType(column[format[key]])
            elif key == "gtf_decorators": # special exceptions for gtf files
                gtf = column[format["gtf_decorators"]].strip()
                for item in gtf.split(";"):
                    if item:
                        item = item.strip()
                        key = item.split(" ")[0]
                        value = item.split(" ")[1].strip("\"")
                        d[key] = self._guessDataType(value)
        return(d)

    def save(self, filename=None, compressed=False):
        """
        **Purpose**

            Save the genelist as a binary representation.
            This is guaranteed to be available for all geneList representations, with
            the only exception being the delayedlists. As that wouldn't
            make any sense as delayedlists are not copied into memory.
            
            You can use this method to cache the file. It's particularly useful for large files
            that get processed once but are then used a lot. 
            
            loading the list back into memory is relatively quick.
            
            list = glload("path/to/filename.glb")
            
            I generally used extension is glb. Although you can use
            whatever you like.

        **Arguments**

            filename
                filename (and path, if you like) to save the file to 

            compressed (Optional, default=False)
                use compression (not currently implemented)

        **Result**

            returns None
            Saves a binary representation of the geneList

        """
        assert filename, "no filename specified"

        oh = open(filename, "wb")
        if compressed:
            config.log.warning("compression not currently implemented, saving anyway")
            cPickle.dump(self, oh, -1)
        else:
            cPickle.dump(self, oh, -1)
        oh.close()
        config.log.info("Saved binary version of list: '%s'" % filename)
