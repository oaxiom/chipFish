"""
track.py tester code.

Part of Chip-Fish

Tests that track is performing accurately.

TODO:
-----
. some of the tests are not very extensive yet.

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../"))

import track
from glbase_wrapper import location

class Test_Track_Function(unittest.TestCase):
    def setUp(self):
        self.t = track.track(filename="test.trk", new=True, name="Test Track")
        self.t.add_location(location(loc="chr1:10-20"))
        self.t.add_location(location(loc="chr1:10-20")) # test duplicates
        self.t.add_location(location(loc="chr1:21-22")) # test 1's
        self.t.add_location(location(loc="chr1:23-23")) # test single outside
        self.t.add_location(location(loc="chr1:1-100")) # test massive span
        self.t.add_location(location(loc="chr1:9-19")) # inside test
        self.t.add_location(location(loc="chr1:15-15")) # 1bp inside
        self.t.add_location(location(loc="chr1:19-25")) # over right border
        self.t.add_location(location(loc="chr1:5-11")) # over left border
        self.t.add_location(location(loc="chrX:1-1000")) # letter chromsome
        self.t.add_location(location(loc="chr2:2-2000")) # other numeric chr
        self.t.add_location(location(loc="chr3:1-2"), strand="-") # test strand
        self.t.add_location(location(loc="chr10:1-2"), strand="-") # test strand
        self.t.add_location(location(loc="chr100:1-2"), strand="-") # test strand

    def test_get_array(self):
        a = self.t.get_array(location(loc="chr1:10-20"))
        self.assertEqual(str(a), "array('i', [5, 5, 4, 4, 4, 5, 4, 4, 4, 5, 4])")

    def test_get_reads(self):
        """
        this wont work after a reload as the strand stuff comes back as a unicode string (u'+')
        """
        a = self.t.get_reads(location(loc="chr1:1-100"))
        self.assertEqual(str(a), "[(10, 20, '+'), (10, 20, '+'), (21, 22, '+'), (23, 23, '+'), (1, 100, '+'), (9, 19, '+'), (15, 15, '+'), (19, 25, '+'), (5, 11, '+')]") # all seqs.
        a = self.t.get_reads(location(loc="chrX:10-20"))
        self.assertEqual(str(a), "[(1, 1000, '+')]")
        a = self.t.get_reads(location(loc="chr3:1-5"))
        self.assertEqual(str(a), "[(1, 2, '-')]")
        a = self.t.get_reads(location(loc="chr1:10-20"))
        self.assertEqual(str(a), "[(10, 20, '+'), (10, 20, '+'), (1, 100, '+'), (9, 19, '+'), (15, 15, '+'), (19, 25, '+'), (5, 11, '+')]")
        a = self.t.get_reads(location(loc="chr1:15-15"))
        self.assertEqual(str(a), "[(10, 20, '+'), (10, 20, '+'), (1, 100, '+'), (9, 19, '+'), (15, 15, '+')]")

    def save_reload(self):
        pass

    def test_read_extend(self):
        a = self.t.get_array(location(loc="chr1:10-20"), read_extend=1)
        self.assertEqual(str(a), "array('i', [5, 5, 5, 4, 4, 5, 5, 4, 4, 5, 5])") # .
        a = self.t.get_array(location(loc="chr1:10-20"), read_extend=2)
        self.assertEqual(str(a), "array('i', [5, 5, 5, 5, 4, 5, 5, 5, 4, 5, 5])") # .

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Track_Function)
    unittest.TextTestRunner(verbosity=2).run(suite)
