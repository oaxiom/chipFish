"""

bookmarks.py

implemented as a sql database.

"""

import sys, os, sqlite3

from glbase_wrapper import location

from userconfig import userconfig

class bookmarks:
    def __init__(self, genome):
        """
        **Purpose**
            Deal with the bookmarks

        **Arguments**
            genome
                genome should be the name of the genome, e.g. "mm8", "Hg18"
        """
        self.__connection = None
        self.__cursor = None

        # get the location of the user data:
        self.userpath = str(userconfig()) # userdata tests for read/write access to this dir.

        # Is the genome dir already available?
        self.genome_path = os.path.join(self.userpath, genome.lower())
        if not os.path.exists(self.genome_path):
            os.mkdir(self.genome_path)

        # see if there is already a bookmark db there:
        self.boomark_db_filename = os.path.join(self.genome_path, "bookmarks.db")
        if not os.path.exists(self.boomark_db_filename):
            self.__setup_database()
        else:
            # I need to test if the db is readable.
            # in case it is mangled in some way.
            self.__connection = sqlite3.connect(self.boomark_db_filename)
            self.__cursor = self.__connection.cursor()
            # test if table 'marks' is availble.
            try:
                self.__cursor.execute("SELECT * FROM marks")
            except:
                self.__setup_database() # not available, so I make the db anyways.

        # open a connection to the db:
        if not self.__connection:
            self.__connection = sqlite3.connect(self.boomark_db_filename)
        if not self.__cursor:
            self.__cursor = self.__connection.cursor()

    def __setup_database(self):
        """
        (Internal)
        setup a vanilla db for the bookmarks
        """
        self.__connection = sqlite3.connect(self.boomark_db_filename)
        self.__connection.text_factory = sqlite3.OptimizedUnicode # support non-latin texts

        if not self.__cursor:
            self.__cursor = self.__connection.cursor()

        self.__cursor.execute("CREATE TABLE marks (i INTEGER PRIMARY KEY AUTOINCREMENT, loc TEXT, chrom TEXT, left INT, right INT, notes TEXT)")

        self.__connection.commit()
        return(None)

    def _getdbpointer(self):
        """
        (Private)
        return a connection to the underlying sqldb
        """
        return(self.__connection)

    def add_bookmark(self, location, notes):
        """
        **Purpose**
            Add a bookmark

        **Arguments**
            location
                a location datatype

            notes
                other notes for the bookmark
                (suggested you use a list of nearby genes)

        **Returns**
            Nothing.
            Remeber to regenerate the bookmark menu
            using self.get_menus()
        """
        # check for a duplicate in the db?
        # DO I already have this entry?
        self.__cursor.execute("SELECT * FROM marks WHERE loc=? AND notes=?",
            (str(location), notes)) # Don't add duplicates.

        if not self.__cursor.fetchall():
            self.__cursor.execute("INSERT INTO marks VALUES (NULL, ?, ?, ?, ?, ?)",
                (str(location), location["chr"], location["left"], location["right"], notes))
            self.__connection.commit()

    def get_all_bookmarks(self):
        """
        **Purpose**
            returns a list of all of the bookmarks in the db.

        **Arguments**
            None

        **Results**
            returns a list of dictionaries of the form
            [{"location": <location>, "notes": <string>} .. n]
        """
        self.__cursor.execute("SELECT * FROM marks")
        # sensibly format:

        return [
            {"loc": location(loc=i[1]), "notes": i[5]}
            for i in self.__cursor.fetchall()
        ]

    def del_bookmark(self, location, notes):
        pass

if __name__ == "__main__":
    # to go into a test suite later.

    b = bookmarks("mm8")
    b.add_bookmark(location(loc="chr17:15064087-15088782"), "Dll1")
    b.add_bookmark(location(loc="chr17:15064087-15088782"), "Dll3")
    print(b._get_all_bookmarks())
