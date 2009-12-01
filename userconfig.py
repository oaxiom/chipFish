"""

userconfig.py

os-specific dealing with paths for local storage of config and data
on a per-user basis.

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
        self.path = "None"
        self.userpath = os.path.expanduser("~")

        if os.name == "posix":
            self.path = self.__posix(self.userpath)
        elif os.name == "nt":
            raise ErrorOSNotSupported
            self.path = self.__win(self.userpath)
        elif os.name == "mac":
            raise ErrorOSNotSupported
            self.path = self.__mac(self.userpath)
        else:
            raise ErrorOSNotSupported

    def __posix(self, userpath):
        # see if a local path already exists:
        config_path = os.path.join(userpath, ".chipfish_userdata")

        if os.path.exists(config_path):
            # see if we can read/write to it.
            if os.access(config_path, os.R_OK | os.W_OK):
                return(config_path) # all seems well.
        else: # dir does not exist
            # see if I own and can write here.
            if os.access(userpath, os.R_OK | os.W_OK):
                try:
                    os.mkdir(config_path)
                except:
                    raise ErrorUserDataReadOnly
                return(config_path)
            else:
                raise ErrorUserDataNotFound

    def __win(self, userpath):
        return("")

    def __mac(self, userpath):
        return("")

    def __str__(self):
        """
        (Override)
        return the path to the per-user data.
        """
        return(self.path)

if __name__ == "__main__":
    # this goes in a test suite later

    a = userconfig()
    print a # should print the
