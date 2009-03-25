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

class bBox:
    """
    A bounding box class for click zones.
    Used in the genome browser view.
    """
    def __init__(self, _dimensions, _location, _object, _type):
        self.dim = _dimensions # a (fourpule)
        self.ptr = _object  # a pointer to some sort of object (or None)
        self.type = _type
        self.loc = _location # location of the item;

    def getObject(self):
        return(self.ptr)

    def getDimensions(self):
        return(self.dim)

    def getTopLeft(self):
        return((self.loc[0], self.loc[1]))

    def getBottomRight(self):
        return((self.loc[0]+self.dim[2], self.loc[1]+self.dim[3]))

    def setLocation(self, _location):
        self.loc = _location

    def setDimensions(self, _dimensions):
        self.dim = _dimensions

    def getRectangle(self):
        """
        returns a rect suitable for pygame.draw.rect()
        """
        print "l:", self.loc
        print "d:", self.dim
        return( (self.loc[0]+self.dim[0], self.loc[1]+self.dim[1],  self.dim[2], self.dim[3]))

    def collideA(self,loc_tuple):
        """
        test for a collision;
        just returns true or false, doens't return the object info.
        use collideB if you want the object's information.
        """
        x = loc_tuple[0]
        y = loc_tuple[1]

        if (x > self.loc[0]+self.dim[0]) and (x < self.loc[0]+self.dim[2]):
            if (y > self.loc[1] + self.dim[1]) and (y < self.loc[1] + self.dim[3]):
                return(True)
        return(False)

    def collideB(self,loc_tuple):
        """
        test for a collision;
        return a dict containing the type and the object
        """
        x = loc_tuple[0]
        y = loc_tuple[1]

        if (x > self.loc[0]+self.dim[0]) and (x < self.loc[0]+self.dim[2]):
            if (y > self.loc[1] + self.dim[1]) and (y < self.loc[1] + self.dim[3]):
                return( {"type": self.type, "object": self.ptr} )
        return(False)
