"""
track pass through for chipFish, see the accompanying readme.txt for details.

(c) 2009-2012 oAxiom

"""

import sys, os, numpy

#sys.path.append("..") # get the parent options
import opt

import glbase # import like this to get around namespace issues

class flat_track(glbase.flat_track):
    """
    chipFish override of vanilla glbase track.py

    In chipFish I want to store some information about the current track
    and
    """
    # descriptions, text for localisation
    __doc__ = "Overriden: Not Present"
    __tooltype__ = "Vanilla track"
    _default_draw_type = "graph"
    _available_draw_types = ("graph", "bar")

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
                strand is an option, but is ignored

            resolution
                the bp resolution to use for the array.

        **Results**
            returns a Numpy array, or a dictionary.
        """
        data = self.get(loc, **kargs)
        
        # flat_tracks do not respect a resolution argument.
        # This should probably be ported back to glbase
        
        newa = numpy.zeros(int(len(data)/resolution))
        
        for i in xrange(len(newa)):
            newa[i] = data[int(i*resolution)]
        
        if type in self._available_draw_types:
            return(newa)

        raise AssertionError, "draw mode not available"

    # gui stuff.
    # available for the gui on this class
    __gui__avail__ = {
        "draw type": _available_draw_types
        }


