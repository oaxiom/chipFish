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
from bookmarks import bookmarks

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
        """
        (Event)

        Executed on the initialisation of the form.
        Main entry point for chipFish.
        """
        # errors should be put into a log;
        #sys.stderr = self # silence errors.
        #sys.stdout = self
        print "chipFish (c) 2009, oAxiom, %s" % (time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())) # log.

        # load the gui
        self.res = xrc.XmlResource('gui_parent.xrc')

        # load the mainFrame
        self.main = self.res.LoadFrame(None, "mainFrame")

        # --------------------------------------------------------------
        # These are user state variables, currently hard coded
        # but later must go freefrom.
        # set up and load the genome
        # (In future this will load the initial state here)
        self.g = glload("data/mm8_refGene.glb")
        self.draw = gDraw.gDraw(self.g) # drawer must have access to the genome

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

        # end of user hard-coded segments.
        # --------------------------------------------------------------

        # menu names:
        # Boomarks root:
        # bookmarks
        # bookmarkThisLocation ID_BOOKMARK_THIS_LOCATION
        # ID_ORGANISE_BOOKMARKS organiseBookmarks
        # findGene
        self.__rebind_bookmarks("mm8") # I have to bodge this as self.g has a mangled name.

        # bind the gPanel;
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.gui_panel = wx.xrc.XRCCTRL(self.main, "gui_panel")
        gDrawPanel = wx.xrc.XRCCTRL(self.main, "gDrawPanel")
        self.draw.bindPanel(gDrawPanel) # bind the drawPanel to gDraw
        draw_panel = self.draw.getPanel()
        sizer.Add(draw_panel, 2, wx.EXPAND, 0) # add the panel to the gui
        gDrawPanel.SetSizer(sizer) # add it to the GUI

        # bind events to the GUI.
        self.Bind(wx.EVT_LEFT_DOWN, self._mouseLeftDown, draw_panel)
        self.Bind(wx.EVT_LEFT_UP, self._mouseLeftUp, draw_panel)
        self.Bind(wx.EVT_BUTTON, self.OnViewLeft, xrc.XRCCTRL(self.main, "butViewLeft"))
        self.Bind(wx.EVT_BUTTON, self.OnViewRight, xrc.XRCCTRL(self.main, "butViewRight"))
        self.Bind(wx.EVT_BUTTON, self.OnViewBigRight, xrc.XRCCTRL(self.main, "butViewBigRight"))
        self.Bind(wx.EVT_BUTTON, self.OnViewBigLeft, xrc.XRCCTRL(self.main, "butViewBigLeft"))
        self.Bind(wx.EVT_BUTTON, self.OnButZoomIn, xrc.XRCCTRL(self.main, "butZoomIn"))
        self.Bind(wx.EVT_BUTTON, self.OnButZoomOut, xrc.XRCCTRL(self.main, "butZoomOut"))
        self.Bind(wx.EVT_BUTTON, self.OnGotoLocationEdit, xrc.XRCCTRL(self.main, "butGoToLoc"))

        # menu events should get by ID: ID_ABOUT
        self.Bind(wx.EVT_MENU, self.OnMenuAbout, id=xrc.XRCID("about"))
        self.Bind(wx.EVT_MENU, self.OnMenuQuit, id=xrc.XRCID("quit"))
        self.Bind(wx.EVT_MENU, self.OnMenuAddBookmark, id=xrc.XRCID("bookmarkThisLocation"))
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
        self.lastValidLocation = self.draw.getLocation()
        if redrawDisplay:
            self.draw.forceRedraw()
        # text boxes.
        self.textGoToLoc.SetValue(str(location(chr=self.draw.chromosome, left=self.draw.lbp, right=self.draw.rbp)))

    def __mouseLeftDown(self, event):
        """
        (Interal)
        Deal with mouse down presses
        """
        pass

    def __mouseLeftUp(self, event):
        """
        (Internal)
        Deal with mouse up presses.
        """
        pass

    def __rebind_bookmarks(self, genome_name):
        """
        (Internal)
        get the bookmarks based on the genome name/identifier
        and bind them into the bookmarks mend
        """
        # get the bookmarks
        self.__bookmarks = bookmarks(genome_name)
        self.__bookmark_lookup = {}

        marks = self.__bookmarks.get_all_bookmarks()

        # the menu to load the bookmarks into
        menu_head = xrc.XRCID("bookmarks")
        menu_head = xrc.XRCCTRL(self.main, "bookmarks")
        menu_head = self.main.GetMenuBar().GetMenu(3) # Bookmarks is 3... IS this stable?

        # clear the previous bookmarks:
        for i in xrange(50000, 50100): # this is a limit to the number of bookmarks...
            # could you even display 100 bookmarks?
            try:
                menu_head.Delete(id=i)
            except:
                break # reached the end.

        if not marks:
            # load a 'None' into the menu
            menu_head.Append(50000, text="None") # how to make inactive?
        else:
            for index, item in enumerate(marks):
                if index > 100:
                    break # I don't support more than 100 at the current time...
                # 50000 + index = ID, well above any other menu ID.
                menu_head.Append(50000+index, text="%s - (%s)" % (item["loc"], item["notes"]))

                # bind a event
                self.Bind(wx.EVT_MENU,
                    self.OnBookmark,
                    id=50000+index)
                self.__bookmark_lookup[50000+index] = item["loc"] # maintain a lookback of location versus Id.

    #------------------------------------------------------------------
    # Events

    def OnBookmark(self, event):
        # get the parent.
        self.draw.setLocation(loc=self.__bookmark_lookup[event.GetId()])
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

    def OnMenuAddBookmark(self, event):
        nearby_features = self.g.getAllDrawableFeaturesInRange(self.draw.getLocation())
        if len(nearby_features) < 8: # truncate the list if it gets too long.
            notes = ", ".join([x["name"] for x in nearby_features])
        else:
            notes = ", ".join([x["name"] for x in nearby_features[0:7]]+["..."])
        self.__bookmarks.add_bookmark(self.draw.getLocation(), notes)
        self.__rebind_bookmarks("mm8") # bodge for now.
# end of mainFrame
