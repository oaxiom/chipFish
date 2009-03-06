"""
chromosome, part of chipFish

(c) 2008-2009 oAxiom 

Not for distribution.

"""

import error, env, gDraw

# these are more 'how to draw' than actual genuine types. 
# for example Chip-seq coords could be Point, whilst
validFeatureTypes = [
	"Point", "Span", "Peak", "Expression" # generic types
	"Gene", "ncRNA", "histone" # specific architecture
	]

def feature:
	def __init__(self, featureFile, featureType, gDraw_ptr):
		"""
		This is my internal implemetation of features
		"""
		self.type = feactureType
		self.file = featureFile
		self.draw = gDraw_ptr
	
	def _exportFeatures(self):
		"""
		Pickle the feature
		"""
		pass
		
	def _importFeatures(self):
		"""
		unPickle the feature
		"""
		pass

	def drawFeature(self, location):
		"""
		draw this feature across location to surface.
		interfaces with gDraw
		"""
		
		# collect data;
		data = 
		
		if self.type = "Point":
			self.draw.drawPoint(surface, location, data)
		elif self.type == "Span":
			self._drawSpan(surface, location, data)
		elif self.type == "Gene":
			self._drawGene(surface, location, data)
