"""
error, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class containter for error handling.

. remove the error class, and raise exceptions.
. assertions to add
"""

import sys, os

# A generic error class.
class error:
    def __init__(self, message, bFatal=False):
        print "Error: ", message
        if bFatal:
            sys.exit()

# ----------------------------------------------------------------------
# Generic Assertion Error:

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

# ---------------------------------------------------------------------
# Exceptions

class ErrorCairoDraw(Exception):
    """
    Error: Attempting to draw to a non-ready Cairo surface.
    Attributes:
        message - an explanation of the error.
    """
    def __init__(self):
        self.message = "Error: Attempted to draw on an uninitialised Cairo Surface."
    def __str__(self):
        return(repr(self.message))

class ErrorCairoAcquireDevice(Exception):
    """
    Error: Attempting to draw to a non-ready Cairo surface.
    Attributes:
        message - an explanation of the error.
    """
    def __init__(self):
        self.message = "Error: Failed to acquire a Cairo Surface."
    def __str__(self):
        return(repr(self.message))

class ErrorInvalidGeneDefinition(Exception):
    """
    Error: Attempting to draw to a non-ready Cairo surface.
    Attributes:
        message - an explanation of the error.
    """
    def __init__(self):
        self.message = "Error: Failed to acquire a Cairo Surface."
    def __str__(self):
        return(repr(self.message))

class ErrorInvalidChromosome(Exception):
    """
    Error: Invalid Chromosome name, only 1-999 and X, Y, M are valid chromsome names.
    Attributes:
        message - an explanation of the error.
    """
    def __init__(self):
        self.message = "Error: Invalid Chromosome name, only 1-999 and X, Y, M are valid chromsome names."
    def __str__(self):
        return(repr(self.message))

class ErrorLibrarySQLite(Exception):
    """
    Error: SQL backend not available.
    """
    def __init__(self):
        self.message = "Error: Database SQLite nota found or not available."
    def __str__(self):
        return(repr(self.message))
