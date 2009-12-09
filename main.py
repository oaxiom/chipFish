"""
main, part of chipFish

(c) 2008-2010 oAxiom

Not for distribution.

"""

import wx

from gui import cfApp

# genome and thread set-up's here;

app = cfApp()
app.MainLoop()

"""
TODO:
-----

. subclass ruler so I can add scale bars
. speed up the bar draws.
"""
