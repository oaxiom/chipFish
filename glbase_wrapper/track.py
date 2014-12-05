"""
track pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os

sys.path.append("..") # get the parent options
import opt
import glbase # import like this to get around namespace issues

class track(glbase.track):
    """
    chipFish override of vanilla glbase track.py

    In chipFish I want to store some information about the current track
    and
    """
    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Vanilla track"
    _default_draw_type = "graph"
    _available_draw_types = ("graph", "graph_split_strand", "bar", "kde_graph")

    def get_data(self, type, loc, strand=None, resolution=1, norm_by_lib_size=False, **kargs):
        """
        **Purpose**
            get data from the track.
            entry point to get the graph
            for track.py it just maps onto get_array()

        **Arguments**
            loc
                location span to get the array for.

            strand (True|False)
                If true return a dictionary {"+": array(), "-": array()}
                type should also == "graph_split_strand"?

            resolution
                the bp resolution to use for the array.

        **Results**
            returns a Numpy array, or a dictionary.
        """
        if type in self._available_draw_types:
            if norm_by_lib_size:
                trk = self.get(loc, resolution=resolution, **kargs) * 100000000 # Quick hack to get it back to ints.
                return(trk / self.get_total_num_reads())
            if type == "kde_graph":
                return(self.get(loc, resolution=resolution, kde_smooth=True, **kargs))
            else:
                if not strand:
                    return(self.get(loc, resolution=resolution, **kargs))
                else:
                    return({"+": self.get(loc, resolution=resolution, strand="+", **kargs),
                        "-": self.get(loc, resolution=resolution, strand="-", **kargs)})
        raise AssertionError, "draw mode not available"

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        "draw type": _available_draw_types
        }


