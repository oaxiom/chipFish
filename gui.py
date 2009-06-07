"""
gui, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for the gui,

TODO:
. introduce a log class

"""

import math, sys, time

import opt, gDraw

#from genome import *
from error import *

import wx
from wx import xrc

# ----------------------------------------------------------------------
# Use a modified start-up of glbase.
# ----------------------------------------------------------------------

# Find glbase
sys.path.append("../../glbase/")

# modify the start-up of glbase
import config # some of the module names will clash...
config.REnv.libraries_to_test = []

# get the libraries I want.
from flags import *
from genelist import genelist, load
from microarray import microarray
from genome import genome
from taglist import taglist
from seqfile import seqfile
from peaklist import peaklist
import utils

# ----------------------------------------------------------------------
# Main GUI
# ----------------------------------------------------------------------

class cfApp(wx.App):
    """
    main frame of the application
    """
    def OnInit(self):
        # errors should be put into a log;
        #sys.stderr = self # silence errors.
        #sys.stdout = self
        print "chipFish (c) 2009, oAxiom, %s" % (time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())) # log.

        # set up and load the genome
        # (In future this will load the initial state here)

        self.g = load("data/mm8.glb")

        # load the gui
        self.res = xrc.XmlResource('gui_parent.xrc')

        # load the mainFrame
        self.main = self.res.LoadFrame(None, "mainFrame")

        # bind the gPanel;
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.draw = gDraw.gDraw(self.g)
        self.gui_panel = wx.xrc.XRCCTRL(self.main, "gui_panel")
        gDrawPanel = wx.xrc.XRCCTRL(self.main, "gDrawPanel")
        self.draw.bindPanel(gDrawPanel) # bind the drawPanel to gDraw
        draw_panel = self.draw.getPanel()
        sizer.Add(draw_panel, 2, wx.EXPAND, 0) # add the panel to the gui
        gDrawPanel.SetSizer(sizer) # add it to the GUI
        
        self.draw.setLocation("3", 153772000, 153850000) # set the location of the genome.
        self.draw.setLocation("3", 153844000, 153850000) # set the location of the genome.

        # bind events to the GUI.
        self.Bind(wx.EVT_LEFT_DOWN, self._mouseLeftDown, draw_panel)
        self.Bind(wx.EVT_LEFT_UP, self._mouseLeftUp, draw_panel)
        self.Bind(wx.EVT_BUTTON, self.OnViewLeft, wx.xrc.XRCCTRL(self.main, "butViewLeft"))
        self.Bind(wx.EVT_BUTTON, self.OnViewRight, wx.xrc.XRCCTRL(self.main, "butViewRight"))
        self.Bind(wx.EVT_BUTTON, self.OnViewBigRight, wx.xrc.XRCCTRL(self.main, "butViewBigRight"))
        self.Bind(wx.EVT_BUTTON, self.OnViewBigLeft, wx.xrc.XRCCTRL(self.main, "butViewBigLeft"))
        self.Bind(wx.EVT_BUTTON, self.OnButZoomIn, wx.xrc.XRCCTRL(self.main, "butZoomIn"))
        self.Bind(wx.EVT_BUTTON, self.OnButZoomOut, wx.xrc.XRCCTRL(self.main, "butZoomOut"))

        # get changable elements from the gui and store them locally.
        # (See _updateDisplay())
        self.textGoToLoc = wx.xrc.XRCCTRL(self.main, "textGoToLoc")
        #print "textGoToLoc", self.textGoToLoc

        self.main.Show()
        self.main.Maximize(True)
        print "End %s" % time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())
        return(True)

    #------------------------------------------------------------------
    # Internal
    def _updateDisplay(self, redrawDisplay=True, event=None):
        """
        (Internal)
        update the gui forms and changeable elements.
        redrawDisplay will not redraw the gui, it just redraws the
        gDrawPanel
        """
        self.lastValidLocation = utils.formatLocation(self.draw.chromosome, self.draw.lbp, self.draw.rbp)
        if redrawDisplay:
            self.draw.forceRedraw()
        # text boxes.
        self.textGoToLoc.SetValue(utils.formatLocation(self.draw.chromosome, self.draw.lbp, self.draw.rbp))

    def _mouseLeftDown(self, event):
        """
        (Interal)
        Deal with mouse down presses
        """
        print "mouse down"
        #if self.D.isColliding():
        #   print "!!!"
        #self.drag = True

    def _mouseLeftUp(self, event):
        """
        (Internal)
        Deal with mouse up presses.
        """
        pass

    #------------------------------------------------------------------
    # Events

    def OnBigViewLeft(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        # move eft depending upon the scale.
        self.draw.setLocation(self.D.chromosome, self.D.lbp - 10000, self.D.rbp - 10000)
        self._updateDisplay()
        event.Skip()

    def OnViewLeft(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self.draw.setLocation(self.D.chromosome, self.D.lbp - 1000, self.D.rbp - 1000)
        self._updateDisplay()
        event.Skip()

    def OnJumpToGenomicLoc(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self._updateDisplay()
        event.Skip()

    def OnViewRight(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self.draw.setLocation(self.D.chromosome, self.D.lbp + 1000, self.D.rbp + 1000)
        self._updateDisplay()
        event.Skip()

    def OnViewBigRight(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self.draw.setLocation(self.D.chromosome, self.D.lbp + 10000, self.D.rbp + 10000)
        self._updateDisplay()
        event.Skip()

    def OnViewBigLeft(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self.draw.setLocation(self.D.chromosome, self.D.lbp - 10000, self.D.rbp - 10000)
        self._updateDisplay()
        event.Skip()

    def OnButZoomIn(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        if self.draw.delta > self.D.w:
            self.draw.setLocation(self.draw.chromosome, self.draw.lbp + 1000, self.draw.rbp - 1000)
        self._updateDisplay()
        event.Skip()

    def OnButZoomOut(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        self.draw.setLocation(self.draw.chromosome, self.draw.lbp - 1000, self.draw.rbp + 1000)
        self._updateDisplay()
        event.Skip()

    def OnGotoLocationEdit(self, event): # wxGlade: frame_mainFrame_parent.<event_handler>
        loc = utils.getLocation(self.textGoToLoc.GetValue())
        if loc:
            self.draw.setLocation(loc["chr"], loc["left"], loc["right"])
        else:
            pass
        self._updateDisplay()
        event.Skip()
# end of mainFrame


if __name__ == "__main__":
    app = cfApp()
    app.MainLoop()
