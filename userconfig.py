"""

userconfig.py

os-specific dealing with paths for local storage of config and data
on a per-user basis.

The loading and saving of items is only partially started.

"""

import sys, os, platform

from error import ErrorUserDataNotFound, ErrorUserDataReadOnly, ErrorOSNotSupported

class userconfig:
    def __init__(self):
        """
        **Purpose**
            generate and deal with paths to find an appropriate
            location to save per-user data.
        """
        self.__userpath = os.path.expanduser("~")

        func_map = {"posix": self.__posix,
            "mac": self.__mac}
            #"win": self.__win, # not currently implemented.

        if os.name not in func_map:
            raise ErrorOSNotSupported

        func_map[os.name]()

        # at this point several paths are now valid:
        # self.__userpath
        # self.__userdatapath = the user data path (on unix ~/.chipfish_userdata)
        # self.__desktop = path to the desktop

        self.__load()

    def __posix(self):
        # see if a local path already exists:
        config_path = os.path.join(self.__userpath, ".chipfish_userdata")
        self.__desktop = os.path.join(self.__userpath, "Desktop") # is this correct for all distributions? okay for Fedora and Ubuntu...

        if os.path.exists(config_path):
            # see if we can read/write to it.
            if os.access(config_path, os.R_OK | os.W_OK):
                self.__userdatapath = config_path # all seems well.
        else: # dir does not exist
            # see if I own and can write here.
            if os.access(self.__userpath, os.R_OK | os.W_OK):
                try:
                    os.mkdir(config_path)
                except:
                    raise ErrorUserDataReadOnly
                self.__userdatapath = config_path
            else:
                raise ErrorUserDataNotFound

    def __win(self):
        """
        I don't have a windows machine.
        """
        return("")

    def __mac(self):
        """
        (Internal)
        As far as I can tell, best practice on MacOSX is identical to
        posix
        """
        self.__posix()

    def __str__(self):
        """
        (Override)
        return the path to the per-user data.
        """
        return(self.__userdatapath)

    def __load(self):
        """
        (Internal)
        Load in the user-data, and if none found set to
        sensible defaults.
        """
        userdata_storage_path = os.path.join(str(self), ".userdata")

        if not os.path.exists(userdata_storage_path):
            # set up sensible defaults.
            self.__data = {"last_png_file": self.__desktop}
        else:
            # load the data in and check for sanity
            pass

    def save(self):
        """
        (Internal)
        save the user data to a file.
        """
        pass

    def __getitem__(self, key):
        """
        (Override)
        Confer
        a = userconfig["last_png_path"]
        """
        return(self.__data[key])

    def __setitem__(self, key, value):
        """
        (Override)
        make properties writable.

        userconfig["last_png_path"] = ""
        """
        self.__data[key] = value

if __name__ == "__main__":
    # this goes in a test suite later

    a = userconfig()
    print a # should print the current userpath.
