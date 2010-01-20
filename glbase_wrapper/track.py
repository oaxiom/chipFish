"""
track pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009 oAxiom

"""

import sys, os

sys.path.append("..") # get the parent options
import opt

sys.path.append(opt.path.glbase_package)
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
    _default_draw_type = "graph_split_strand"
    _available_draw_types = ("graph", "graph_split_strand", "bar")

    def get_data(self, type, loc, strand=None, resolution=1, **kargs):
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

            resolution
                the bp resolution to use for the array.

        **Results**
            returns a Numpy array, or a dictionary.
        """
        if type in self._available_draw_types:
            if not strand:
                return(self.get_array(loc, resolution=resolution, **kargs))
            else:
                return({"+": self.get_array(loc, resolution=resolution, strand="+"),
                    "-": self.get_array(loc, resolution=resolution, strand="-")})
        raise AssertionError, "draw mode not available"

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        "draw type": _available_draw_types
        }


