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
        Problems all over here: Why does the start-up not match other start ups?
        bad bad bad

        **Arguments**

            None

        """

        format_override = default
        genelist.__init__(self)
        if kargs.has_key("format"):
             format_override = kargs["format"]

        # defaults:

        for k in kargs:
            if k == "name":
                self.name = kargs["name"]
            elif k == "filename":
                self.loadCSV(filename=kargs["filename"], format=format_override)

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
            if item.has_key("loc"):
                loc = deepcopy(item["loc"])
                loc.pointify()
                new_item = deepcopy(item)
                new_item["loc"] = loc
                newl.linearData.append(new_item)
        newl._optimiseData()
        self._history.append("Pointified (locations converted to a single base pair, in the middle of the coordinates)")
        return(newl)

    def frequencyAgainstArray(self, **kargs):
        """
        Draw an peaklist and compare against an array.
        Draws a three panel

        **Arguments**

            microarray

                must be a microarray-like object containing a "conditions" key, and other
                condition descriptors.

            filename

                a file name to the iamge as, including the path.

        **Result**

            returns True if succesful or False if unsuccesful.

            plots a multi-panel image of the array and a moving average plot of the density of
            chip-seq peaks. Saves the image to filename.
        """
        # reqd args enforcer.
        req_args = ["filename", "microarray"]
        res = [i for i in req_args if i not in kargs]
        if res:
            print "Error: in method frequencyAgainstArray(), these required arguments are missing: %s" % (", ".join(res))
            sys.quit()
        # end of the enforcer.

        kargs["array_data"] = kargs["microarray"]
        kargs["peak_data"] = self
        kargs["row_names"] = kargs["microarray"]["name"]
        kargs["col_names"] = kargs["microarray"].getConditionNames()
        if not kargs.has_key("bracket"): kargs["bracket"] = [0, 1]

        actual_filename = self.draw._heatmap_and_plot(**kargs)

        if not config.SILENT: print "Info: Saved heatmap and moving average plot to %s" % actual_filename
        return(True)
