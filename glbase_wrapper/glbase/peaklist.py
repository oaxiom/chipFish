"""
An inherited version of geneList that has extra methods for
dealing with chip-seq peaks.
"""

import config, utils
from genelist import genelist
from draw import draw
from copy import deepcopy
from flags import *
from errors import AssertionError

class peaklist(genelist):
    def __init__(self, **kargs):
        """
        **Purpose**

        Problems all over here: Why does the start-up not match other start ups?
        bad bad bad

        **Arguments**
            filename (Optional, if not specified, returns an empty peaklist)
                filename of the file to load.

            format (Optional, defaults to sniffer)
                An optional format specifier

            Name (Optional, defaults to filename)
                override the name, if not specified, use the file name instead.

        **Returns**
            a valid peaklist.
        """

        format_override = default # from flags.
        genelist.__init__(self)
        if "format" in kargs and kargs["format"]:
             format_override = kargs["format"]

        # defaults:

        if "filename" in kargs and kargs["filename"]:
            self.loadCSV(filename=kargs["filename"], format=format_override)

        if "name" in kargs and kargs["name"]:
            self.name = kargs["name"]

    def __repr__(self):
        return("glbase.peaklist")

    def pointify(self):
        """
        Convert all of the loc coordinates to a single point, centred
        around the middle of the coordinates

        Your list must have a 'loc' key

        **Arguments**

            None

        **Result**

            Returns a new list with 'pointified' coordinates - the coordinates
            are now a single base pair centred around the middle of the
            coordinates.
        """
        newl = self.__copy__()
        newl.linearData = []
        for item in self:
            if "loc" in item:
                loc = deepcopy(item["loc"])
                loc.pointify()
                new_item = deepcopy(item)
                new_item["loc"] = loc
                newl.linearData.append(new_item)
        newl._optimiseData()
        self._history.append("Pointified (locations converted to a single base pair, in the middle of the coordinates)")
        return(newl)

    def frequencyAgainstArray(self, filename=None, microarray=None, **kargs):
        """
        Draw an peaklist and compare against an array.
        Draws a three panel figure showing an array heatmap, the binding events
        and a moving average of the binding events.

        **Arguments**

            microarray

                must be a microarray-like object containing a "conditions" key, and other
                condition descriptors.

            filename

                a file name for the image including the path.

            window (Optional)

                size of the sliding window to use.

            bracket (Optional)

                bracket the data within a range (e.g. (0,2))

        **Result**

            returns True if succesful or False if unsuccesful.

            plots a multi-panel image of the array and a moving average plot of the density of
            chip-seq peaks. Saves the image to filename.
        """
        valid_args = ["filename", "microarray", "window", "bracket"]
        for key in kargs:
            assert key in valid_args, "unrecognised argument %s" % key

        assert filename, "must specify a filename to save as"
        assert microarray, "must provide a microarray list"

        # reload/override not really a good way to do this...
        # I should reload a new dict... As I may inadvertantly override another argument?
        kargs["filename"] = filename
        kargs["array_data"] = microarray
        kargs["peak_data"] = self
        kargs["row_names"] = microarray["name"]
        kargs["col_names"] = microarray.getConditionNames()
        if not "bracket" in kargs: kargs["bracket"] = [0, 1]

        actual_filename = self.draw._heatmap_and_plot(**kargs)

        if not config.SILENT: print "Info: Saved heatmap and moving average plot to %s" % actual_filename
        return(True)
