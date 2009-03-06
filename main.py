"""
main, part of chipFish

(c) 2008-2009 oAxiom 

Not for distribution.

"""

import wx

from gui import mainFrame

# genome and thread set-up's here;

chipFishApp = wx.PySimpleApp(0)
wx.InitAllImageHandlers()
mainFrame = mainFrame(None, -1, "")
chipFishApp.SetTopWindow(mainFrame)
mainFrame.Show()
chipFishApp.MainLoop()
