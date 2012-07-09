
"""

This should be how to use the chipFish from the command line.

And output svgs etc.

"""


from app import app
from glbase_wrapper import location

a = app()
a.startup("../Tracks/stat3_tracks.txt")

a.draw.setLocation(loc=location(loc="chr11:100747467-100813886"))
a.draw.exportImage("Stat3.svg", type="svg")

a.draw.setLocation(loc=location(loc="chr11:100747467-100813886"))
a.draw.exportImage("Stat3.svg", type="svg")

a.draw.setLocation(loc=location(loc="chr11:100747467-100813886"))
a.draw.exportImage("Stat3.svg", type="svg")

a.draw.setLocation(loc=location(loc="chr11:100747467-100813886"))
a.draw.exportImage("Stat3.svg", type="svg")