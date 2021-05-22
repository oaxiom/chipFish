"""
main, part of chipFish

(c) 2008-2010 oAxiom

Not for distribution.

"""


import wx, profile

import opt, sys
from gui import cfApp

app = cfApp()
# genome and thread set-up's here;

if not opt.debug.profile:
    app.MainLoop()
else:
    import cProfile, pstats
    cProfile.run("app.MainLoop()", "profile.pro")
    sys.stdout = open("log.log", "w")
    p =  pstats.Stats("profile.pro")
    p.strip_dirs().sort_stats("time").print_stats()

"""
TODO:
-----

. subclass ruler so I can add scale bars
. speed up the bar draws.
"""
