"""
An inherited version of geneList that has extra methods for
dealing with chip-seq peaks.
"""

import config, utils
from genelist import genelist
from draw import draw
from copy import deepcopy
from flags import *
from errors import AssertionError, NoMatchingKeysError, BadOperationError

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

            name (Optional, defaults to filename)
                override the name, if not specified, use the file name instead.

        **Returns**
            a valid peaklist.
        """
        valid_args = ["filename", "format", "name"]
        for key in kargs:
            assert key in valid_args, "unrecognised argument %s" % key

        format_override = default # from flags.
        genelist.__init__(self)
        if "format" in kargs and kargs["format"]:
             format_override = kargs["format"]

        # defaults:

        if "filename" in kargs and kargs["filename"]:
            self.loadCSV(filename=kargs["filename"], format=format_override)
            config.log.info("Loaded '%s' found %s items" % (kargs["filename"], len(self.linearData)))

        if "name" in kargs and kargs["name"]:
            self.name = kargs["name"]

        config.log.debug("Initialised peaklist '%s'" % self.name)

    def __repr__(self):
        return("glbase.peaklist")

    def pointify(self, key="loc"):
        """
        Convert all of the loc coordinates to a single point, centred
        around the middle of the coordinates

        Uses a 'loc' key as the default location to pointify

        **Arguments**

            key (default = "loc")
                specify a location key to pointify.
                defaults to loc


        **Result**

            Returns a new list with 'pointified' coordinates - the coordinates
            are now a single base pair centred around the middle of the
            coordinates.
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.__copy__()
        for item in newl:
            loc = item["loc"].pointify()
        newl._optimiseData()
        config.log.info("pointified peaklist %s" % self.name)
        self._history.append("Pointified (locations converted to a single base pair, in the middle of the coordinates)")
        return(newl)

    def frequencyAgainstArray(self, filename=None, microarray=None, **kargs):
        """
        Draw a peaklist and compare against an array.
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
                Defaults to 10% of the length of the microarray.
                If you specify a number here it is in items on the array
                and not a percentage.

            bracket (Optional)

                bracket the data within a range (e.g. (0,2))

            match_key (Optional)
                the key to match between the array and the peaklist.
                Goes throgh a cascade until it finds one of:
                ["refseq" > "entrez" > "loc" > "tss_loc" >
                "array_systematic_name" > "name"]

            tag_key (Optional)
                use a 'tag score' or some other key:value pair to use
                as a metric for the peaklist.

                By default the list is trated as a binary list, either
                it is or is not a binding site. More sophisticated analysis
                may use the number of seq reads, or fold change over
                background or some other metric to give a more analogue-like
                score. You can specify the score to use using this argument

        **Result**

            returns True if succesful or False if unsuccesful.

            plots a multi-panel image of the array and a moving average plot of the density of
            chip-seq peaks. Saves the image to filename.
        """
        valid_args = ["filename", "microarray", "window", "bracket", "match_key", "tag_key"]
        for key in kargs:
            assert key in valid_args, "unrecognised argument %s" % key

        assert filename, "must specify a filename to save as"
        assert microarray, "must provide a microarray list"

        # work out a matching key or just use the user-supplied one.
        match_key = None
        if "match_key" in kargs and kargs["match_key"]:
            match_key = kargs["match_key"]
        else:
            keys1 = microarray.getKeys()
            keys2 = self.getKeys()
            keys_to_try = ["refseq", "entrez", "loc", "tss_loc", "array_systematic_name", "name"]
            for k in keys_to_try:
                if k in keys1 and k in keys2:
                    match_key = k
                    break

        if not match_key:
            raise NoMatchingKeysError, ("peaklist.frequencyAgainstArray()", None)

        tag_key = None
        if "tag_key" in kargs and kargs["tag_key"]:
            tag_key = kargs["tag_key"]

        # reload/override not really a good way to do this...
        # I should reload a new dict... As I may inadvertantly override another argument?
        kargs["match_key"] = match_key
        kargs["filename"] = filename
        kargs["array_data"] = microarray
        kargs["peak_data"] = self
        kargs["row_names"] = microarray["name"]
        kargs["col_names"] = microarray.getConditionNames()
        kargs["use_tag_score"] = tag_key
        if not "bracket" in kargs: kargs["bracket"] = [0, 1]

        actual_filename = self.draw._heatmap_and_plot(**kargs)

        config.log.info("Saved heatmap and moving average plot to %s" % actual_filename)
        return(True)

