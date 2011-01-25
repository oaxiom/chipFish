"""
config.py

config must be imported before any other glbase library.

"""

# config cannot import any other glbase module

import os, logging

# container for global environment variables.

# The current version of GLbase
VERSION = "0.166.hg"
DATE = "17th Jan 2011"

SILENT = False # set this to True to silence all glbase output.
DEBUG = True

# flags for the availability of three core libraries. [Deprecated now?]
MATPLOTLIB_AVAIL = False
NUMPY_AVAIL = False
SCIPY_AVAIL = False

# Some simple options for printing genelists
NUM_ITEMS_TO_PRINT = 3 # number of items to print by default.
PRINT_LAST_ITEM = True

DEFAULT_DPI = 150
DEFAULT_DRAWER = "png" # slated for deprecation.

# set up the logger here.
# this needs to be moved to log.py
# You can access it using config.log()
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-8s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=os.path.join(os.path.expanduser("~"), "glbase.log"),
                    filemode='w')
console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(logging.Formatter('%(levelname)-8s: %(message)s'))# # console
logging.getLogger('').addHandler(console)

# use config.log. ... () to get to the logger
log = logging.getLogger('glbase')

# helpers [Should be deprecated. Call config.log.<level>() to call info
info = log.info
warning = log.warning
debug = log.debug
error = log.error

if SILENT:
    log.setLevel(logging.CRITICAL) # not acutally silenced...

def set_log_level(level):
    """
    Change the logging level for the on-screen logger.
    the console logger is not affected by this call.
    """
    level_map = {None: logging.CRITICAL, # not correct?
        "info": logging.INFO,
        "debug": logging.DEBUG}

    assert level in level_map, "no valid level used to set the logger, valid modes are 'info' or 'debug'"

    log.setLevel(level_map[level])
    return(True)
