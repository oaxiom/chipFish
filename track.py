"""
track, part of chipFish

(c) 2008-2009 oAxiom 

Not for distribution.

Class container for the track data

"""

import utils, error

class track:
	def __init__(self, name):
		self.name = name
		
	def setTag(self, tag, value):
		"""
		meta data about the track
		"""
		pass

	def setData(self):
		"""
		I don't unneccesarily load the data until it is required
		"""
		pass
	
	def _loadData(self):
		"""
		genuinely load the data
		"""
		pass
	
	def _unloadData(self):
		"""
		unload the data as not required.
		"""
		pass
		
	def unpackData(self, location):
		loc = utils.getLocation(location)
		if not loc:
			error.error("Location not valid", False)
			return(False)
			
			 
