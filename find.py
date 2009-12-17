"""
find.py
"""

import wx
from wx import xrc

showable_keys = ("name", "entrez", "refseq", "loc")

class findDialog:
    def __init__(self, parent_frame):
        # load the frame
        self.res = xrc.XmlResource('gui_parent.xrc')
        self.main = self.res.LoadDialog(None, "CustomFindDialog")
        self.parent_frame = parent_frame

        # get relevant elements:
        self.__m_textFindText = wx.xrc.XRCCTRL(self.main, "m_textFindText")
        self.__m_gauge = wx.xrc.XRCCTRL(self.main, "m_gaugeSearch")
        self.__m_listBox = wx.xrc.XRCCTRL(self.main, "m_listBox")
        self.__m_butGoTo = wx.xrc.XRCCTRL(self.main, "m_butGoTo")

        # bind events.
        self.main.Bind(wx.EVT_BUTTON, self.OnStartSearch, xrc.XRCCTRL(self.main, "m_butSearchStart"))
        self.main.Bind(wx.EVT_BUTTON, self.OnGoTo, xrc.XRCCTRL(self.main, "m_butGoTo"))

        self.main.Show()

    def OnStartSearch(self,event):
        """
        Start searching for the text in parent_frame.g
        """
        search_string = self.__m_textFindText.GetValue()

        results = self.parent_frame.g.find(search_string, self.__m_gauge)

        # load the results into the ListCtrl

        if not results:
            self.__m_listBox.Set([])
            return(None)

        # look at the first item to get useful keys:
        mykeys =[]
        for key in results[0]:
            if key in showable_keys:
                mykeys.append(key)

        newl = []
        for item in results:
            newl.append(", ".join([str(item[key]) for key in mykeys])) # reformat

        self.__m_listBox.Set(newl)
        self.__last_found = results

    def OnGoTo(self, event):
        """
        move the genomic view to the location specified by location
        """
        # get current item;
        citem = self.__m_listBox.GetSelections()

        if len(citem):
            citem = self.__last_found[citem[0]]

        self.parent_frame._updateDisplay(loc=citem["loc"]) # eh... probably shouldn't reach into its guts like that...
