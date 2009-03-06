"""
gui, part of chipFish

(c) 2008-2009 oAxiom 

Not for distribution.

Class container for constants and other data

"""

import wx, guiWxParent

class mainFrame(guiWxParent.frame_mainFrame_parent):    
    def __init__(self, *args, **kwds): 
        guiWxParent.frame_mainFrame_parent.__init__(self, *args, **kwds) # inherit
        
        self.dchandle = wx.PaintDC
        
        self.Maximize(True)

