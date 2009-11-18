"""
gui, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for the gui,

TODO:
. introduce a log class

"""

import math, sys, time

import opt, gDraw, glbase_wrapper

from error import *

import wx
from wx import xrc

# ----------------------------------------------------------------------
# Use a modified start-up of glbase.
# ----------------------------------------------------------------------

# Find glbase
sys.path.append(opt.path.glbase_wrapper)
from glbase_wrapper import *

from track import track

# ----------------------------------------------------------------------
# Main GUI
# ----------------------------------------------------------------------

class cfApp(wx.App):
    def OnInit(self):
        # errors should be put into a log;
        #sys.stderr = self # silence errors.
        #sys.stdout = self
        print "chipFish (c) 2009, oAxiom, %s" % (time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())) # log.

        # set up and load the genome
        # (In future this will load the initial state here)
        self.g = load("data/mm8_refGene.glb")

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

        self.draw.bindTrack(track(filename="data/NSMash1_new.trk", name="NS5 Mash1 ChIP-seq"))
        self.draw.bindTrack(track(filename="data/SpMash1_new.trk", name="Spinal Cord Mash1 ChIP-seq"))
        self.draw.bindTrack(track(filename="data/TcMash1_new.trk", name="Telencephalon Mash1 ChIP-seq"))
        self.draw.bindTrack(track(filename="data/NS_H3K4me3.trk", name="NS5 H3K4me3"))
        #self.draw.bindTrack(track(filename="data/NS_H3K27me3.trk", name="NS5 H3K27me3"))
        #self.draw.bindTrack(track(filename="data/NS_H3K36me3.trk", name="NS5 H3K36me3"))
        #self.draw.bindTrack(track(filename="data/ES_H3K36me3.trk", name="ES H3K4me3"))
        #self.draw.bindTrack(track(filename="data/MEF_H3K36me3.trk", name="MEF H3K4me3"))

        #self.draw.setLocation("6", 122666976, 122685608) # Nice view of Nanog
        #self.draw.setLocation("1", 3001251, 3001551) # testing the ChIP-seq track
        self.draw.setLocation("17", 15074087, 15084782) # Interesting view of the ChIP-seq (Dll1?) chr17:15,064,087-15,088,782
        #self.draw.setLocation("7", 28010095, 28012095) # Dll3
        #self.draw.setLocation("3", 34850525, 34853525) # Sox2
        #self.draw.setLocation("16", 91152965, 91156965) # Olig1
        #self.draw.setLocation("5", 140855715, 140873715) # Lnfg

        # bind events to the GUI.
        self.Bind(wx.EVT_LEFT_DOWN, self._mouseLeftDown, draw_panel)
        self.Bind(wx.EVT_LEFT_UP, self._mouseLeftUp, draw_panel)
        self.Bind(wx.EVT_BUTTON, self.OnViewLeft, wx.xrc.XRCCTRL(self.main, "butViewLeft"))
        self.Bind(wx.EVT_BUTTON, self.OnViewRight, wx.xrc.XRCCTRL(self.main, "butViewRight"))
        self.Bind(wx.EVT_BUTTON, self.OnViewBigRight, wx.xrc.XRCCTRL(self.main, "butViewBigRight"))
        self.Bind(wx.EVT_BUTTON, self.OnViewBigLeft, wx.xrc.XRCCTRL(self.main, "butViewBigLeft"))
        self.Bind(wx.EVT_BUTTON, self.OnButZoomIn, wx.xrc.XRCCTRL(self.main, "butZoomIn"))
        self.Bind(wx.EVT_BUTTON, self.OnButZoomOut, wx.xrc.XRCCTRL(self.main, "butZoomOut"))
        self.Bind(wx.EVT_BUTTON, self.OnGotoLocationEdit, wx.xrc.XRCCTRL(self.main, "butGoToLoc"))

        # menu events should get by ID: ID_ABOUT
        self.Bind(wx.EVT_MENU, self.OnMenuAbout, id=xrc.XRCID("about"))
        self.Bind(wx.EVT_MENU, self.OnMenuQuit, id=xrc.XRCID("quit"))
        # get changable elements from the gui and store them locally.
        # (See _updateDisplay())
        self.textGoToLoc = wx.xrc.XRCCTRL(self.main, "textGoToLoc")
        #print "textGoToLoc", self.textGoToLoc

        self.main.Show()
        self.main.Maximize(True)
        # force a redraw:
        self._updateDisplay()
        print "End %s" % time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
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
        self.lastValidLocation = location(chr=self.draw.chromosome, left=self.draw.lbp, right=self.draw.rbp)
        if redrawDisplay:
            self.draw.forceRedraw()
        # text boxes.
        self.textGoToLoc.SetValue(str(location(chr=self.draw.chromosome, left=self.draw.lbp, right=self.draw.rbp)))

    def _mouseLeftDown(self, event):
        """
        (Interal)
        Deal with mouse down presses
        """
        pass

    def _mouseLeftUp(self, event):
        """
        (Internal)
        Deal with mouse up presses.
        """
        pass

    #------------------------------------------------------------------
    # Events

    def OnJumpToGenomicLoc(self, event):
        string_loc = self.textGoToLoc.text
        print text_loc
        #self.draw.setLocation(self.
        self._updateDisplay()
        event.Skip()

    def OnViewRight(self, event):
        self.draw.move("right", 10)
        self._updateDisplay()
        event.Skip()

    def OnViewBigRight(self, event):
        self.draw.move("right", 20)
        self._updateDisplay()
        event.Skip()

    def OnViewLeft(self, event):
        self.draw.move("left", 10)
        self._updateDisplay()
        event.Skip()

    def OnViewBigLeft(self, event):
        self.draw.move("left", 20)
        self._updateDisplay()
        event.Skip()

    def OnButZoomIn(self, event):
        self.draw.move("zoomin", 10)
        self._updateDisplay()
        event.Skip()

    def OnButZoomOut(self, event):
        self.draw.move("zoomout", 10)
        self._updateDisplay()
        event.Skip()

    def OnGotoLocationEdit(self, event):
        try:
            loc = location(loc=self.textGoToLoc.GetValue())

            if loc:
                self.draw.setLocation(loc["chr"], loc["left"], loc["right"])
            else:
                raise TypeError
        except (TypeError, IndexError):
            # colour the textBox Red;
            pass
        self._updateDisplay()
        event.Skip()

    # ------------------------------------------------------------------
    # menu events
    def OnMenuAbout(self, event):
        info = wx.AboutDialogInfo()
        info.SetName("ChipFish")
        info.SetDescription("Genomic Data and Analysis, in a salty and vinegary package.\nWrapped in newspaper with mayonnaise")
        info.SetCopyright("(c) 2009-2010 oAxiom")
        info.SetVersion("\n\nchipfish: %s\n glbase: %s" %(opt.generic.VERSION, glbase_wrapper.VERSION))
        info.AddDeveloper("Andrew Hutchins")
        info.SetWebSite("http://www.oaxiom.com")
        wx.AboutBox(info)

    def OnMenuQuit(self, event):
        # check to see if we want to save our work

        # any other user cleanups.

        # and die:
        sys.quit()
# end of mainFrame
