"""
track, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

"""

import cPickle, sys, os, struct, ConfigParser # renamed to configparser in >2.6

from glbase_wrapper import location as location
from glbase_wrapper import positive_strand_labels, negative_strand_labels

from array import array

bin_len = {"i": struct.calcsize("i"), "l": struct.calcsize("l")}

class track:
    """
    track definition, used for things like sequence reads across the genome
    """
    def __init__(self, stranded=True, delayed=True, dir_name="", bin_format="i"):
        """
        **Arguments**
            stranded (True|False)
                store one or both strands?

            delayed (True|False)
                set as a delayed edit?
                (Do not load into memory)

            dir_name (string)
                directory location of the track file.
                only respected if dir_name is set.

        """
        self.delayed = delayed
        self.track = {}
        self.bin_format = bin_format
        if self.delayed:
            self._load(dir_name)

    def add_location(self, loc, strand="+"):
        assert not self.delayed, "the list must not be 'delayed' to constuct the initial track sets"

        # make sure chromosome is available:
        if not self.track.has_key(loc["chr"]):
            self.track[loc["chr"]] = {"+": array(self.bin_format , []), "-": array(self.bin_format , [])}

        # check the stranded nature: make sure we're using the same symbol.
        # (chipFish prefers +/-
        if strand in positive_strand_labels:
            strand = "+"
        elif strand in negative_strand_labels:
            strand = "-"

        # make sure the postion is available
        if len(self.track[loc["chr"]][strand]) < loc["right"]:
            # nope, track is too short
            size_required = loc["right"] - len(self.track[loc["chr"]][strand])
            self.track[loc["chr"]][strand].extend([0 for x in xrange(size_required)])

        for l in xrange(loc["left"], loc["right"]):
            self.track[loc["chr"]][strand][l] += 1

    def __getitem__(self, loc, strand="+"):
        """
        t[location(loc="chr1:10-20")]
        how to do strands though?

        ** slow **
            It may be slow, but is required
        """
        if self.delayed:
            #try:
            self.track[loc["chr"]].seek(loc["left"] * bin_len[self.bin_format])
            readlen = loc["right"] - loc["left"]
            value = self.track[loc["chr"]].read(readlen * bin_len[self.bin_format])
            real = [value[x:x+bin_len[self.bin_format]] for x in xrange(0, readlen, bin_len[self.bin_format])]
            new_a = array(self.bin_format, [])
            for item in real:
                print item
                print struct.unpack(self.bin_format, item)[0]
                print "%x" % struct.unpack(self.bin_format, item)[0]
                new_a.append(struct.unpack(self.bin_format, item)[0])
            return(new_a)
            #except TypeError: # probably a str
            #    print "FF"

        elif not self.delayed:
            # the data is immediatly available.
            try:
                return(self.track[loc["chr"]][strand][loc["left"]:loc["right"]])
            except TypeError: # treat loc as if it was a string:
                loc = location(loc=loc)
                return(self.track[loc["chr"]][strand][loc["left"]:loc["right"]]) # try again...

    def save(self, dir_name):
        """
        save the track as a delayed format set of track files.
        """
        # this is a human readable descriptor file describing the track
        if not os.path.exists(dir_name): # this is dangerous?
            os.mkdir(dir_name)

        descriptor_file = ConfigParser.RawConfigParser()
        descriptor_file.add_section("File Manifest")
        descriptor_file.add_section("Element Sizes")
        for chrom in self.track:
            descriptor_file.set("File Manifest", "name", "%s.trk" % chrom)
            for strand in ["+", "-"]:
                descriptor_file.set("Element Sizes", "%s(%s)" % (chrom, strand), str(len(self.track[chrom][strand])))

        descriptor_file.add_section("Format")
        descriptor_file.set("Format", "binary descriptor", self.bin_format)
        descriptor_file.set("Format", "chromosomes", str([key for key in self.track]))

        oh = open(os.path.join(dir_name, "descriptor.txt"), "w")
        descriptor_file.write(oh)
        oh.close()

        # py2.6 version:
        #with open(os.path.join(dir_name, "descriptor.txt"), "w") as file:
        #    descriptor_file.write(file)

        for chrom in self.track:
            oh = open(os.path.join(dir_name, "%s.trk" % chrom), "wb")

            strand = ["+", "-"]
            for s in strand:
                oh.write("%s\n" % s)
                for b in self.track[chrom]["+"]:
                    oh.write(struct.pack(self.bin_format , b)) # not floatable... # need to have some sort of format tag for the track
                oh.write("\n")
            oh.close()
        return(True)

    def _load(self, dir_name):
        """
        opposite of load, sets up the track so that track[] works okay.
        call from the init, don't call directly.
        """
        # get the file descriptor.
        self.delayed = True

        descriptor_file = ConfigParser.RawConfigParser()
        descriptor_file.read(os.path.join(dir_name, "descriptor.txt"))

        self.bin_format = descriptor_file.get("Format", "binary descriptor")
        chroms = eval(descriptor_file.get("Format", "chromosomes"))

        # some fo the config files are not fully used...
        for chrom in chroms:
            chrom_handle = open(os.path.join(dir_name, "%s.trk" % chrom), "r")
            self.track[chrom] = chrom_handle

if __name__ == "__main__":
    # a few testers to see how it works.

    t = track(delayed=False)
    t.add_location(location(loc="chr1:10-20"))
    t.add_location(location(loc="chr1:11-21"))
    t.add_location(location(loc="chr1:12-22"))
    t.add_location(location(loc="chr1:25-26"))
    t.add_location(location(loc="chr1:25-26"))
    t.add_location(location(loc="chr1:25-26"))
    # really boost one single location:
    for n in xrange(1000):
        t.add_location(location(loc="chr1:11-12"))
        t.add_location(location(loc="chr1:12-13"))

    t.add_location(location(loc="chr1:10-20"), strand="-")
    t.add_location(location(loc="chr1:11-21"), strand="-")
    t.add_location(location(loc="chr1:12-22"), strand="-")
    t.add_location(location(loc="chr1:25-26"), strand="-")

    print t.track
    print t[location(loc="chr1:10-20")]
    print t["chr1:10-20"]
    t.save("data/test/")

    #print t[1] # error

    t = track(delayed=True, dir_name="data/test/")
    print t[location(loc="chr1:10-20")]
    print t["chr1:10-20"]
