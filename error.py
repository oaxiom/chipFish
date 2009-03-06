"""
error, part of chipFish

(c) 2008-2009 oAxiom 

Not for distribution.

Class container for utilities for conversion etc...

"""

import sys, os

class error:
	def __init__(self, message, bFatal=False):
		print "Error: ", message
		if bFatal:
			sys.exit()		
		
