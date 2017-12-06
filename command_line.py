
"""

This should be how to use the chipFish from the command line.

And output svgs etc.

"""

import os
from app import app
from glbase_wrapper import location

a = app()
a.startup(os.path.expanduser("~/Projects/AP2aGCM/trks/track_list.txt"))
oh = open(os.path.expanduser("~/Projects/AP2aGCM/trks/For raw signal.txt"), "rU")

for lin in oh:
    tt = lin.strip().split("\t")
    
    print(tt)
    if tt[0] != "gname":
    
        a.draw.setLocation(loc=location(loc=tt[1]).expand(10000))
        a.draw.exportImage(os.path.expanduser("~/Projects/AP2aGCM/trks/%s_%s.png" % (tt[0], str(tt[1].replace(":", "-")))), type="png")