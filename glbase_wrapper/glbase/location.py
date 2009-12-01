"""

location.py

part of glbase.

This class is an internal class that implemnts a more convenient way to manipulate
genomic coordiantes.

TODO:
. add a 'in' code clause e.g.:
    if 1000 in location: (see if 1000 > left & < right)
    if a_location in b_location: (exectute a collide())

"""

import copy

from errors import AssertionError

class location:
    def __init__(self, loc=None, chr=None, left=None, right=None):
        if isinstance(loc, location):
            # It's actually already a loc.
            # I want to copy it and leave.
            self._loc_string = "chr%s:%s-%s" % (loc["chr"].strip("chr").lower(), loc["left"], loc["right"])
            self.loc = copy.copy(loc.loc)
        else:
            if loc:
                self._loc_string = loc.lower()
                t = self._loc_string.split(":")
                self.loc = {"chr": t[0].strip("chr").upper(), "left":int(t[1].split("-")[0]), "right":int(t[1].split("-")[1])}
            else:
                #self._loc_string = "chr%s:%s-%s" % (chr.lower().strip("chr"), left, right)
                #t = self._loc_string.split(":")
                self.loc = {"chr": chr.strip("chr").upper(), "left": int(left), "right": int(right)}

    def __repr__(self):
        self._loc_string = self._merge(self.loc) # only update when accessed.
        return("<location %s>" % (self._loc_string))

    def __len__(self):
        self._loc_string = self._merge(self.loc) # only update when accessed.
        return(len(self._loc_string)) # this should be length of span!?!?!

    def split(self, value=None):
        # ignores the 'value' argument completely and returns a three-ple
        return( (self.loc["chr"], self.loc["left"], self.loc["right"]) )

    def _merge(self, loc_dict):
        try:
            return("chr%s:%s-%s" % (loc_dict["chr"].strip("chr"), loc_dict["left"], loc_dict["right"]))
        except: # chr possibly sets of strings ... etc.
            return("chr%s:%s-%s" % (loc_dict["chr"], loc_dict["left"], loc_dict["right"]))
        return(False) # bugged out;

    def __getitem__(self, key):
        if key == "string":
            self._loc_string = self._merge(self.loc) # only update when accessed.
            return(self._loc_string)
        if key == "dict":
            return(self.loc)
        return(self.loc[key])

    def __setitem__(self, key, value):
        self.loc[key] = value

    def __str__(self):
        self._loc_string = self._merge(self.loc) # only update when accessed.
        return(self._loc_string)

    """
    these methods below should copy the location and send a modified version back.
    """
    def expand(self, base_pairs):
        self.loc["left"] -= base_pairs
        self.loc["right"] += base_pairs
        return(self)

    def expandLeft(self, base_pairs):
        self.loc["left"] -= base_pairs
        return(self)

    def expandRight(self, base_pairs):
        self.loc["right"] += base_pairs
        return(self)

    def shrink(self, base_pairs):
        self.loc["left"] += base_pairs
        self.loc["right"] -= base_pairs
        return(self)

    def shrinkLeft(self, base_pairs):
        self.loc["left"] += base_pairs
        return(self)

    def shrinkRight(self, base_pairs):
        self.loc["right"] -= base_pairs
        return(self)

    def pointify(self):
        centre = (self.loc["left"] + self.loc["right"]) / 2
        self.loc = {"chr": self.loc["chr"], "left": centre, "right": centre}
        return(self)

