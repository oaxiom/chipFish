"""
db, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

interface to the db backend.

abstracted to add in/out other db backends.
based on SQL and using SQLite for this particular implmenetation.

At the moment, it's all done with a bunch of adhoc csv's :(

"""

import sqlite

class db:
    def __init__(self, name, database, **kargs):
        """
        valid kargs
        readonly=True
        """
        for item in kargs:
            print item

        self.name = name
        self.database = database

    def getQuery(self, sql_string):
        pass

    def setQuery(self, sql_string):
        pass

    def commitChanges(self):
        pass

    def rollBackChanges(self):
        pass

