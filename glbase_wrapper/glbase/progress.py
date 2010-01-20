
"""

a simple progress bar for some of the longer functions expected to take a long-time.

Inspired from the progress bar in urlgrabber.

"""

from __future__ import division

import sys

import config

class progressbar:
    def __init__(self, maximum_value, output=sys.stderr):
        """
        (Override)
        Initialise the progress bar

        **Arguments**

            maximum_value (Required)
                the maximum_value to move the bar up to.

            output (Optional, defaults to stderr)
                the output device to use.
        """
        self.maximum = maximum_value
        self.__writer = output
        self.__barwidth = 30 # bar_width in characters. This may need to change on later computers
                             # with larger terminals
        self.__last_percent = 0 # only print if __last_pecent is incremented.

    def update(self, new_value):
        """
        Update progress meter with new_value

        **Arguments**

            new_value (Required)
                should be some number between 0 .. maximum_value
        """
        percent_done = int(((new_value+1) / self.maximum) *100)

        if self.__last_percent < percent_done:
            t_percent_done = int(((new_value+1) / self.maximum) * self.__barwidth)

            bar = "".join(["=" for x in xrange(t_percent_done)] + ["-" for x in xrange(self.__barwidth-t_percent_done)])
            if not config.SILENT: self.__writer.write("\r[%s] %s%% (%s/%s)" % (bar, percent_done, new_value, self.maximum))
            self.__last_percent = percent_done

        if new_value+1 >= self.maximum: # if the last line, reset the console so the result overprints the progress bar.
            if not config.SILENT: self.__writer.write("\r") # pad out to overwrite the previous bar.
            if not config.SILENT: self.__writer.write("\r                                                 ") # pad out to overwrite the previous bar.
            if not config.SILENT: self.__writer.write("\r") # pad out to overwrite the previous bar.

if __name__ == "__main__":
    """
    A quick tester.
    """
    import time

    p = progressbar(100)
    for i in xrange(100):
        time.sleep(0.1)
        p.update(i)