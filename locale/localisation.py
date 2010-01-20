"""

locale.py

part of chipFish

handles and loads the localisation.txt file.

chipFish then looks up (by string) in what is essentially a huge dictionary
of language data.

"""

import sys, os

from wx import xrc

supported_languages = frozenset(["en-gb"])

class locale:
    """
    **Purpose**
        emulate a dictionary
    """
    def __init__(self, language):
        assert language in supported_languages, "language: %s not currently supported" % language

        print "Info: Register Localisation: %s" % language

        self.locale = language
        self.__path_to_language_file = os.path.join(sys.path[0], "locale", self.locale, "language.txt")
        self.__data = {}
        self.__load()

        print "Info: Registered %s translatable items" % len(self.__data)

    def __load(self):
        assert os.path.exists(self.__path_to_language_file), "language file not found"

        oh = open(self.__path_to_language_file, "rU")
        for line in oh:
            line = line.strip("\n").strip("\r")
            if len(line) and line[0] not in ["#", "", " ", "\n"]:
                head = line.split("\t")[0]
                tail = line.split("\t")[1]

                self.__data[head] = tail
        oh.close()

    def __getitem__(self, key):
        return(self.__data[key])

    def load_main_gui(self, gui_handle):
        """
        **Purpose**
            The main gui frame of glbase contains a lot of locale
            settings. Instead of stuffing all the locale stuff into there,
            instead I stuff it here. It's going to have to make a mess
            somewhere. Adn here is better as locale is likely to
            be relatively lightweight.

        **Arguments**
            gui_handle
                the handle to the main chipFish gui from which I can
                call xrc.XRCCTRL() on.
        """
        elementsToModify = ["butZoomIn", "butZoomOut", "butGoToLoc",
            "butChrUp", "butChrDown"]

        for element in elementsToModify:
            xrc.XRCCTRL(gui_handle, element).SetLabel(self.__data[element])

        return(True)

        # the menus are a little different:
        elementsToModify = ["file", "about", "help"]

        menu_head = gui_handle.GetMenuBar()

        for element in elementsToModify:
            item = menu_head.FindItem(element)
            print item
            if item:
                menu_item = menu_head.GetMenu(item)
                menu_item.SetLabel(element)

    def load_search_gui(self):
        """
        **Purpose**
            see load_main_gui() for the justification.
        """
        pass
