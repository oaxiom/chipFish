"""
errors.py

clean-up code and helpers for catching and dealing with errors.

"""

import csv

from data import typical_headers, ignorekeys
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
        print "Error: %s" % (message)
        if not config.DEBUG: # this only works in Py3k right?
            self.__traceback__ = None

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
        print
        print "Error: Unrecognised CSV"
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
        print "Format Specifier: %s" % (" ".join(["%s:%s" % (key, format[key]) for key in format]))
        print "Expected Format, based on the format specifier:"
        oh.close()

        # This is a safe version of loadCSV() that intelligently fails.

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
                                        if not d.has_key(key):
                                            d[key] = {}
                                        if isinstance(format[key], dict) and (format[key].has_key("code")):
                                            # a code block insertion goes here - any valid lib and one line python code fragment
                                            # store it as a dict with te key "code"
                                            d[key] = eval(format[key]["code"])
                                        else:
                                            # attempt to coerce things as an int.
                                            # no need to do the cooersion here, just turn it into a str()
                                            d[key] = str(column[format[key]])
                                    except:
                                        d[key] = "mangled"
                            print "%s" % (" ".join(["%s:%s" % (key, d[key]) for key in d]))
                            if index > 3:
                                break
        else:
            print "  No specified format (glbase will guess)"

        print "-----------------------"
        print "Error: %s" % (message)
        print

class LibraryNotFoundError(Exception):
    """
    Error
        A required library is not found.
    """
    def __init__(self, message, traceback):
        print "Error: %s" % (message)

class UnrecognisedFileFormatError(Exception):
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
        print
        print "Error: Unrecognised file format"
        print "-----------------------"
        print "Diagnostic:"
        print "0:", oh.readline().rstrip("\r\n")
        print "1:", oh.readline().rstrip("\r\n")
        print "2:", oh.readline().rstrip("\r\n")
        print "3:", oh.readline().rstrip("\r\n")
        print "Format Specifier: %s" % (" ".join(["%s:%s" % (key, format[key]) for key in format]))
        print "-----------------------"
        print "Error: %s" % (message)
        print
