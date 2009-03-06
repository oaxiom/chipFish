"""
gDraw, part of chipFish

(c) 2008-2009 oAxiom 

Not for distribution.

NOTES:
------
. gDraw is a wrapper for whatever vector-draw backend is in use.
. The vector back-end is not normally available outside of gDraw as its 
	implementation is likely to change.
. DXF/PS exporter?
. I really want to use Cairo here;

"""

import wx.lib.wxcairo

class gDraw:
    def __init__(self, surface):
        pass

	def setViewPort(self, maxx, maxy): 
		"""
		set the size of the viewport;
		"""
		self.maxx = maxx
		self.maxy = maxy
		self.aspect = abs(maxx / maxy)
	
	def setScale(self, scale):
		"""
		scalar 1 - 100
		"""
		self.scale = scale
	
	def setDrawAttribute(self, attribute, value):
		"""
		setDrawAttribute
		"""
		pass

    def drawChr(self, location):
        """
        draw the basic chromosome 
        """

	def drawPoint(self, location, data):
		"""
		drawFeatures of the type "Point"
		data should be a list of coordinates within the location span
		"""
		pass
	
	def drawSpan(self, location, data):
		"""
		drawFeatures of the type "Span"
		data should be a list of coordinates of the form chrX:left-right within the location span
		"""
		pass
	
	def drawGene(self, location, data):
		"""
		drawFeatures of the type "Gene"
		data should be a list of the form:
		[chr, genespan_left, genespan_right, exon_coords(local), strand, cds_start, cds_end]
		"""
		pass
		
	def exportPostScript(self, _filename):
		pass
		
	def exportDXF(self, _filename):
		pass
		
	def exportPNG(self, _filename):
		pass
	
	def exportBMP(self, _filename):
		pass
		
	def exportJPG(self, _filename):
		pass	
