
import copy, cPickle

import config
from errors import AssertionError, UnRecognisedCSVFormatError, UnrecognisedFileFormatError, ArgumentError

class _base_genelist:
    def __init__(self):
        """
        (Internal)
        This is the base derived class for all genelists.
        It contains methods available to all implementations of genelist.
        """
        self.name = None
        self.linearData = None

    def __repr__(self):
        return("<base genelist class>")

    def __nonzero__(self):
        """
        Fixes:
        if genelist_object:
            Now True
        """
        return(True)

    def __copy__(self):
        """
        (Override)
        Confer copy to mean a deep copy as opposed to a shallow copy.
        """
        return(copy.deepcopy(self))

    def __shallowcopy__(self):
        """
        (New)

        Some weird behaviour here, I know, this is so I can still get access to
        the shallow copy mechanism even though 90% of the operations are copies.
        """
        return(copy.copy(self))

    def save(self, filename=None, **kargs):
        """
        **Purpose**

            Save the genelist as a binary representation.
            This is guaranteed to be available for all geneList representations, with
            the only exception being the delayedlists. As that wouldn't
            make any sense as delayedlists are not copied into memory
            this method is used in caching the file.
            use list = glload("path/to/filename.glb")
            The generally used extension is glb. Although you can use
            whatever you like.

        **Arguments**

            filename
                path to the file (including path)

            compressed
                use compression (not currently implemented)

        **Result**

            returns None
            Saves a binary representation of the geneList

        """
        valid_args = ["filename", "compressed"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError, (self.save, k)

        assert filename, "no filename specified"

        compressed = False
        if "compressed" in kargs:
            compressed = kargs["compressed"]
        if compressed:
            config.log.warning("compression not currently implemented, saving anyway")

        oh = open(filename, "wb")
        if compressed:
            cPickle.dump(self, oh, -1)
        else:
            cPickle.dump(self, oh, -1)
        oh.close()
        config.log.info("Saved Binary version of list: '%s'" % filename)
