
"""

This should be how to use the chipFish from the command line.

And output svgs etc.

"""

import os
from app import app
from glbase_wrapper import location

a = app()
a.startup(os.path.expanduser("~/Desktop/Projects/coSTAT3/5.ChIP-seq/stat3_trks/stat3_tracks.txt"))

a.draw.setLocation(loc=location(loc="chr11:100747467-100813886"))
a.draw.exportImage("Stat3.svg", type="svg")

a.draw.setLocation(loc=location(loc="chr11:117801020-117857062"))
a.draw.exportImage("Socs3.svg", type="svg")

a.draw.setLocation(loc=location(loc="chr7:20376837-20425079"))
a.draw.exportImage("Bcl3.svg", type="svg")

a.draw.setLocation(loc=location(loc="chr2:167734770-167805825"))
a.draw.exportImage("Ptpn1.svg", type="svg")