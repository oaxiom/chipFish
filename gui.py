"""
gui, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for the gui, inherits from guiWxParent, auto generated
from wxGlade

"""

import wx, guiWxParent, math, sys, gDraw

class mainFrame(guiWxParent.frame_mainFrame_parent):
    def __init__(self, *args, **kwds):
        """
        The mainFrame for the app.
        """
        #sys.stderr = self # silence errors.
        #sys.stdout = self
        guiWxParent.frame_mainFrame_parent.__init__(self, None, -1) # inherit
        self.Maximize(True)

        # set-up the gDraw Cairo panel
        sizer = wx.BoxSizer(wx.VERTICAL)
        self._Drawer = gDraw.gDraw()
        self._Drawer.bindPanel(self.gDrawPanel)
        sizer.Add(self._Drawer.getPanel(), 2, wx.EXPAND, 0)
        self.gDrawPanel.SetSizer(sizer)
        self.Fit()

if __name__ == "__main__":
    app = wx.App()
    f = mainFrame(None, -1, "")
    f.Show()
    app.MainLoop()
