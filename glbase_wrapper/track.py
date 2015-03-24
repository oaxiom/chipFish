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

    def get_data(self, type, loc, strand=None, resolution=1, norm_by_lib_size=False, norm_factor=1.0, **kargs):
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
        if type not in self._available_draw_types:
            raise AssertionError, "draw mode not available"

        if type == "kde_graph":
            return(self.get(loc, resolution=resolution, kde_smooth=True, **kargs))
        elif type == 'graph_split_strand':
            return({"+": self.get(loc, resolution=resolution, strand="+", **kargs),
                "-": self.get(loc, resolution=resolution, strand="-", **kargs)})
        else:
            trk = self.get(loc, resolution=resolution, **kargs) 
        
        if norm_factor:
            trk *= norm_factor
            
        if norm_by_lib_size:
            trk /= (self.get_total_num_reads() / 100000000.0)# Quick hack to get it back to ints.
            
        return(trk)

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        "draw type": _available_draw_types
        }


