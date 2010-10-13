"""

location.py

part of glbase.

This class is an internal class that implements a more convenient way to manipulate
genomic coordiantes.

TODO:
. add a 'in' code clause e.g.:
    if 1000 in location: (see if 1000 > left & < right)
    if a_location in b_location: (exectute a collide())

"""

import copy

class location:
    def __init__(self, loc=None, chr=None, left=None, right=None):
        if isinstance(loc, location):
            # It's actually already a loc.
            # I want to copy it and leave.
            self._loc_string = "chr%s:%s-%s" % (loc["chr"].strip("chr").lower(), loc["left"], loc["right"])
            self.loc = copy.copy(loc.loc)
        else:
            if loc:
                self._loc_string = loc.lower().replace(",", "") # ucsc includes commas, remove them so you can cut and paste
                t = self._loc_string.split(":")
                self.loc = {"chr": t[0].strip("chr").upper(), "left":int(t[1].split("-")[0]), "right":int(t[1].split("-")[1])}
            else:
                #self._loc_string = "chr%s:%s-%s" % (chr.lower().strip("chr"), left, right)
                #t = self._loc_string.split(":")
                self.loc = {"chr": str(chr).strip("chr").upper(), "left": int(left), "right": int(right)}

    def __eq__(self, other):
        if other:
            if isinstance(other, str):
                return(str(self) == str(other)) # use string comparison.

            # use a faster ? dict comparison, or throw an exception, as this item probably not a <location>
            if self.loc["chr"] == other.loc["chr"]:
                if self.loc["left"] == other.loc["left"]:
                    if self.loc["right"] == other.loc["right"]:
                        return(True)
        return(False)

    def __nonzero__(self):
        return(True)

    def __repr__(self):
        self._loc_string = self._merge(self.loc) # only update when accessed.
        return("<location %s>" % (self._loc_string))

    def __len__(self):
        self._loc_string = self._merge(self.loc) # only update when accessed.
        # work out the span.
        return(self.loc["right"] - self.loc["left"])

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
        new = copy.deepcopy(self)
        new.loc["left"] -= base_pairs
        new.loc["right"] += base_pairs
        return(new)

    def expandLeft(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["left"] -= base_pairs
        return(new)

    def expandRight(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["right"] += base_pairs
        return(new)

    def shrink(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["left"] += base_pairs
        new.loc["right"] -= base_pairs
        return(new)

    def shrinkLeft(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["left"] += base_pairs
        return(new)

    def shrinkRight(self, base_pairs):
        new = copy.deepcopy(self)
        new.loc["right"] -= base_pairs
        return(new)

    def pointify(self):
        new = copy.deepcopy(self)
        centre = (self.loc["left"] + self.loc["right"]) / 2
        new.loc = {"chr": self.loc["chr"], "left": centre, "right": centre}
        return(new)

    def qcollide(self, loc):
        """
        **Purpose**
            perform a collision with another location object.

        **Returns**
            True or False
        """
        if loc["chr"] != self["chr"]:
            return(False)

        # quickest rejections first;
        if self["right"] < loc["left"]:
            return(False)
        if self["left"] > loc["right"]:
            return(False)

        if self["left"] <= loc["right"] and self["right"] >= loc["right"]:
            return(True) # Bright point is within A, collision

        if self["right"] >= loc["left"] and self["left"] <= loc["left"]:
            return(True) # Bleft point is within A, collision.

        if loc["left"] <= self["right"] and loc["right"] >= self["right"]:
            return(True) # Aright point is within B, collision

        if loc["right"] >= self["left"] and loc["left"] <= self["left"]:
            return(True) # Aleft point is within B, collision.

        return(False)

    def distance(self, loc):
        """
        **Purpose**
            calculate the distance between two locations.

        **Returns**
            an integer indicating the distance, note that
            the chromosomes should be the same or it will raise an
            exception
        """
        assert self["chr"] == loc["chr"], "chromosomes are not the same, %s vs %s" % (self, loc)
        return(self.qdistance(loc))

    def qdistance(self, loc):
        """
        (Internal)
        ignore the assert.
        """
        centreA = (self.loc["left"] + self.loc["right"]) / 2
        centreB = (loc["left"] + loc["right"]) / 2
        return(centreA - centreB)

    def offset(self, base_pairs):
        """
        get a new location offset from the 5' end by n base pairs
        returns a point location.
        """
        new = copy.deepcopy(self)
        new.loc["left"] += base_pairs
        new.loc["right"] = new.loc["left"]
        return(new)
