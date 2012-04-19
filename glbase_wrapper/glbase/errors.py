"""
errors.py

clean-up code and helpers for catching and dealing with errors.

"""

import csv

from data import typical_headers, ignorekeys
from location import location
import config

class AssertionError(Exception):
    """
    Error
        An assertion or requirement for a particular method
        failed. This usually means some sort of category required
        for the method is missing.
    """
    def __init__(self, message):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("%s" % (message))
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

class DelayedListError(Exception):
    """
    Error
        Some sort of unsupported function in a delayedlist.
    """
    def __init__(self, message):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("%s" % (message))
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

class ArgumentError(Exception):
    """
    Error
        An assertion or requirement for a particular method
        failed. This usually means some sort of category required
        for the method is missing.
    """
    def __init__(self, function, argument):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("Function '%s' - argument '%s' not supported" % (function, argument))
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

def guessDataType(value):
    """
    (Internal)

    This is a copy of genelist._guessDataTpye()
    I have to copy it here, as due to some circular imports I can't import
    a genelist

    Take a guess at the most reasonable datatype to store value as.
    returns the resulting data type based on a list of logical cooercions
    (explaines as I fail each cooercion).
    Used internally in _loadCSV()
    I expect this will get larger and larger with new datatypes, so it's here as
    as separate proc.

    Datatype cooercion preference:
    float # this is wrong? try to make a float if has '.' otherwise int.
    int
    location
    string
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

class UnRecognisedCSVFormatError(Exception):
    """
    Error
        The csv is not recognised, and produces an error somewhere inside
        _loadCSV(). Print a selection of output based on what glbase
        expects the CSV to look like, hopefully this explains what may be
        going wrong with the CSV.
    """
    def __init__(self, message, file_handle, format):
        """
        Format and ouput a series of messages so that I can see why the csv is not loading.
        """
        oh = open(file_handle, "rU")
        config.log.error("csv/tsv file did not pass the csv parser")
        print "-----------------------"
        print "CSV Diagnostic:"
        if format.has_key("skiplines"): # skip the lines.
            if format["skiplines"] != -1:
                for n in xrange(format["skiplines"]):
                    oh.readline().rstrip("\r\n")

        print "0:", oh.readline().rstrip("\r\n")
        print "1:", oh.readline().rstrip("\r\n")
        print "2:", oh.readline().rstrip("\r\n")
        print "3:", oh.readline().rstrip("\r\n")
        print "-----------------------"
        print "Format Specifier: %s" % (" ".join(["%s:%s\t" % (key, format[key]) for key in format]))
        print "Expected Format, based on the format specifier:"
        oh.close()

        # This is a safe-ish version of loadCSV() that intelligently fails.

        if not format.has_key("sniffer"):
            oh = open(file_handle, "rU")
            if format.has_key("dialect"):
                reader = csv.reader(oh, dialect=format["dialect"])
            else:
                reader = csv.reader(oh)

            try:
                if format.has_key("skiplines"):
                    skiplines = format["skiplines"]
                else:
                    skiplines = 0 # skip any header row by default.
            except:
                print "Error: End of File" # premature end of file, skip out.
                print "-----------------------"
                print "Error: %s" % (message)
                return

            for index, column in enumerate(reader): # This is cryptically called column, when it is actually row.
                if index > skiplines:
                    if column: # list is empty, so omit.
                        if (not (column[0] in typical_headers)):
                            d = {}
                            for key in format:
                                if not (key in ignorekeys): # ignore these tags
                                    try:
                                        if not key in d:
                                            d[key] = {}
                                        if isinstance(format[key], dict) and "code" in format[key]:
                                            # a code block insertion goes here - any valid lib and one line python code fragment
                                            # store it as a dict with the key "code"
                                            d[key] = eval(format[key]["code"]) # this always fails for some reason...
                                        else:
                                            d[key] = str(column[format[key]])
                                    except:
                                        d[key] = "mangled"
                            print "%s" % (" ".join(["%s:%s" % (key, d[key]) for key in d]))
                            if index > 3:
                                break
        else:
            print "  No specified format (glbase will guess)"

        print "-----------------------"
        config.log.error("End of error output")

class LibraryNotFoundError(Exception):
    """
    Error
        A required library is not found.
    """
    def __init__(self, message):
        config.log.critical("Library Not Found: %s" % (message))

class UnrecognisedFileFormatError(Exception):
    """
    Error
        File is not recognised, but not in loadCSV.
        Just print some diagnostic stuff, but not so fancy as
        UnRecognisedCSVFormatError
    """
    def __init__(self, message, file_handle, format):
        """
        Format and ouput a series of messages so that I can see why the csv is not loading.
        """
        oh = open(file_handle, "rU")
        config.log.critical("Unrecognised file format")
        print "-----------------------"
        print "Diagnostic:"
        print "0:", oh.readline().rstrip("\r\n")
        print "1:", oh.readline().rstrip("\r\n")
        print "2:", oh.readline().rstrip("\r\n")
        print "3:", oh.readline().rstrip("\r\n")
        if "sniffer" in format:
            print "Format Specifier: Sniffer (guess the file format)"
        else:
            print "Format Specifier: %s" % (" ".join(["%s:%s" % (key, format[key]) for key in format]))
        print "-----------------------"
        print "Error: %s" % (message)
        print

class NotSupportedError(Exception):
    """
    Error
        Some Command somewhere is not supported.
        For example, delayedlist.collide() only
        implements "and" and the "logic" commands are not supported.
    """
    def __init__(self, message):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("Not Supported: %s" % (message))
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

class NotImplementedError(Exception):
    """
    Error
        Some things I'm too lazy to implement properly, so I raise an error here.
    """
    def __init__(self, message):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("Not Supported: %s" % (message))
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

class NoMatchingKeysError(Exception):
    """
    Error
        No matching keys available between the two lists.
    """
    def __init__(self, function, argument):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        print "Error: No suitable matching key between the two lists"
        print "       valid keys are: refseq entrez loc tss_loc array_systematic_name"
        print "       You can also specify your own key using the 'match_key' argument"
        print "       Both the microarray and the peaklist must have both keys"
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

class BadOperationError(Exception):
    """
    Error
        act() recieved a bad operation value.
    """
    def __init__(self, function, argument):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        print "Error: A bad operation was sent to act(), or act() was unable"
        print "       to complete the action."
        print "       Bad code: %s" % argument
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

class GlglobDuplicateNameError(Exception):
    """
    Error
        glglob recived two names that are identical.
    """
    def __init__(self, function, argument):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("glglob has two lists with the same name")
        config.log.critical("for glglob to work you must provide lists")
        config.log.critical("with unique names")
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

class FailedToMakeNewDBError(Exception):
    """
    Error
        Failed to make a new track database
    """
    def __init__(self, filename):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("Failed to make the new database")
        config.log.critical("This path is not valid or is inaccesible?")
        config.log.critical("Tried: '%s'" % message)
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None
            
class BadBinaryFileFormatError(Exception):
    """
    Error
        The binary file is not a cPickle
    """
    def __init__(self, filename):
        """
        Output the error message and tidy up the traceback, and perform other stuff.
        """
        config.log.critical("File '%s' is not a glbase binary file" % filename)
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None