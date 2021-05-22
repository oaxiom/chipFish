"""
--------------------------
boundbox, part of chipFish

(c) 2008-2009 oAxiom

This source code is
Not for distribution.
--------------------------

bounding box class for dealing with clicking on the
genome display.

NOTES:
------
. This is a self contained class. Only call the methods.
. Borrowed from 'Conglomerate' of all places!
"""

import constants

class bbox:
    """
    A bounding box class for click zones.
    Used in the genome browser view.
    """
    def __init__(self, location, object, type):
        assert len(location) == 4, "wrong sized dimensions supplied to bBox"
        self.dim = location # a (fourpule)
        self.ptr = object  # a pointer to some sort of object (or None)
        self.type = type

    def getObject(self):
        return(self.ptr)

    def get_dimensions(self):
        return(self.dim)

    def collide(self, loc_tuple):
        """
        test for a collision;
        return a dict containing the type and the object
        """
        x = loc_tuple[0]
        y = loc_tuple[1]

        if (
            (x > self.loc[0] + self.dim[0])
            and (x < self.loc[0] + self.dim[2])
            and (y > self.loc[1] + self.dim[1])
            and (y < self.loc[1] + self.dim[3])
        ):
            return( {"type": self.type, "object": self.ptr} )
        return(False)
