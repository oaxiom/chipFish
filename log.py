"""
Logging module for chipfish.

Basically from glbase.

The standard way to call is:
log.info()
log.critical()
log.error

"""

SILENT= False # silence all

# set up the logger here.
import logging, sys, os

# set up logging to file - see previous section for more details
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s: %(message)s',
                    datefmt='%m-%d %H:%M'),

# use config.log. ... () to get to the logger
log = logging.getLogger('chipFish')

# helpers - a bunch of aliases
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

    assert level in level_map, "no valid level used to set the logger"

    log.setLevel(level_map[level])
    return(True)
