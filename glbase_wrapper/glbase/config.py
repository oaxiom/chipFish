"""
config.py

config must be imported before any other glbase library.

"""

# config cannot import any other module

# container for global environment variables.

# The current version of GLbase
VERSION = "0.141"
DATE = "17th Nov 2009"

# Set me to True to get more output from the scripts.
VERBOSE = False

# Set me to False to silence the output from glbase (not including error messages)
SILENT = False

# set me to true to enable some debug output
DEBUG = False

# flags for the availability of three core libraries.
MATPLOTLIB_AVAIL = False
NUMPY_AVAIL = False
SCIPY_AVAIL = False

NUM_ITEMS_TO_PRINT = 3 # number of items to print by default.
PRINT_LAST_ITEM = True

DEFAULT_DPI = 150
DEFAULT_DRAWER = "png"
