"""
gui, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for the gui, inherits from guiWxParent, auto generated
from wxGlade

TODO:
. introduce a log class

"""

import wx, guiWxParent, math, sys, gDraw, utils, wx.xrc

#logic by itself in module

import wx
from wx import xrc


class MyApp(wx.App):
    """
    main frame of the application. Overrides the guiWxParent
    derived from wxGlade
    """
    def OnInit(self):
        # errors should be put into a log;
        #sys.stderr = self # silence errors.
        #sys.stdout = self
        print "chipFish (c) 2009, oAxiom" # log.
        self.res = xrc.XmlResource('gui_parent.xrc')
        # load the mainFrame
        self.main = self.res.LoadFrame(None, "mainFrame")
        # I have to bind all the events here;


        # bind the gPanel;
        sizer = wx.BoxSizer(wx.VERTICAL)
        self._Drawer = gDraw.gDraw()
        self.gui_panel = id=wx.xrc.XRCCTRL(self.main, "gui_panel")
        print self.gui_panel
        #self.frame.Bind(wx.EVT_BUTTON, self.OnSubmit, id=xrc.XRCID('button'))
        gDrawPanel = wx.xrc.XRCCTRL(self.main, "gDrawPanel")
        self._Drawer.bindPanel(gDrawPanel) # bind the drawPanel to gDraw
        draw_panel = self._Drawer.getPanel()
        sizer.Add(draw_panel, 2, wx.EXPAND, 0) # add the panel to the gui
        gDrawPanel.SetSizer(sizer) # add it to the GUI
        self._Drawer.setLocation(3, 9500, 15000)

        print wx.xrc.XRCCTRL(self.gui_panel, "butViewLeft")

        # bind events to the GUI.
        self.Bind(wx.EVT_LEFT_DOWN, self._mouseLeftDown, draw_panel)
        self.Bind(wx.EVT_BUTTON, self.OnViewLeft, wx.xrc.XRCCTRL(self.gui_panel, "butViewLeft"))
        self.Bind(wx.EVT_BUTTON, self.OnViewRight, wx.xrc.XRCCTRL(self.main, "butViewRight"))
        self.Bind(wx.EVT_BUTTON, self.OnViewBigRight, wx.xrc.XRCCTRL(self.main, "butViewBigRight"))


        # get changable elements from the gui
        self.textGoToLoc = wx.xrc.XRCCTRL(self.main, "textGoToLoc")
        print "textGoToLoc", self.textGoToLoc

        self.main.Show()
        self.main.Maximize(True)
        return(True)

        # set-up the gDraw panel
        #sizer = wx.BoxSizer(wx.VERTICAL)
        #self._Drawer = gDraw.gDraw()
        #self._Drawer.bindPanel(self.gDrawPanel) # bind the drawPanel to gDraw
        #draw_panel = self._Drawer.getPanel()
        #sizer.Add(draw_panel, 2, wx.EXPAND, 0) # add the panel to the gui
        #self.gDrawPanel.SetSizer(sizer)

        # Bind mouse events to gDrawPanel
        #self.gDrawPanel.Connect(wx.EVT_LEFT_DOWN, self._mouseLeftDown)
        #self.Bind(wx.EVT_LEFT_DOWN, self._mouseLeftUp)

        #self.Bind(wx.EVT_LEFT_UP, self._mouseLeftUp, sizer)
        #self.Bind(wx.EVT_LEFT_DOWN, self._mouseLeftDown, draw_panel)

        #self.Fit()

        # set-up other parts of the display.
        #self.textGoToLoc.SetValue(utils.formatLocation(3, self._Drawer.lbp, self._Drawer.rbp))

    #------------------------------------------------------------------
    # Internal
    def _updateDisplay(self, redrawDisplay=True, event=None):
        """
        (Internal)
        update the gui forms
        redrawDisplay will not redraw the gui, it just redraws the
        gDrawPanel
        """
        if redrawDisplay:
            self._Drawer.forceRedraw()
        # text boxes.
        self.textGoToLoc.SetValue(utils.formatLocation(3, self._Drawer.lbp, self._Drawer.rbp))

    def _mouseLeftDown(self, event):
        """
        (Interal)
        Deal with mouse down presses
        """
        print "mouse down"
        event.skip()

    def _mouseLeftUp(self):
        """
        (Internal)
        Deal with mouse up presses.
        """
        pass

    #------------------------------------------------------------------
    # Events

    def OnBigViewLeft(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        # move eft depending upon the scale.
        self._Drawer.setLocation(3, self._Drawer.lbp - 10000, self._Drawer.rbp - 10000)
        self._updateDisplay()
        event.Skip()

    def OnViewLeft(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self._Drawer.setLocation(3, self._Drawer.lbp - 1000, self._Drawer.rbp - 1000)
        self._updateDisplay()
        event.Skip()

    def OnJumpToGenomicLoc(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self._updateDisplay()
        event.Skip()

    def OnViewRight(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self._Drawer.setLocation(3, self._Drawer.lbp + 1000, self._Drawer.rbp + 1000)
        self._updateDisplay()
        event.Skip()

    def OnViewBigRight(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self._Drawer.setLocation(3, self._Drawer.lbp + 10000, self._Drawer.rbp + 10000)
        self._updateDisplay()
        event.Skip()

    def OnButZoomIn(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self._Drawer.setLocation(3, self._Drawer.lbp + 1000, self._Drawer.rbp - 1000)
        self._updateDisplay()
        event.Skip()

    def OnButZoomOut(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self._Drawer.setLocation(3, self._Drawer.lbp - 1000, self._Drawer.rbp + 1000)
        self._updateDisplay()
        event.Skip()

    def OnGotoLocationEdit(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        loc = utils.getLocation(self.textGoToLoc.GetValue())
        if loc:
            self._Drawer.setLocation(3, loc["left"], loc["right"])
        else:
            pass
            # change the format of the text.
        event.Skip()
# end of mainFrame


if __name__ == "__main__":
    app = MyApp()
    app.MainLoop()
