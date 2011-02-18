"""
some choice utils from glbase

"""

import math

def qcollide(Aleft, Aright, Bleft, Bright):
    """
    optimised for speed.
    """
    # quickest rejections first;
    if Aright < Bleft:
        return(False)
    if Aleft > Bright:
        return(False)

    if Aleft <= Bright and Aright >= Bright:
        return(True) # Bright point is within A, collision

    if Aright >= Bleft and Aleft <= Bleft:
        return(True) # Bleft point is within A, collision.

    if Bleft <= Aright and Bright >= Aright:
        return(True) # Aright point is within B, collision

    if Bright >= Aleft and Bleft <= Aleft:
        return(True) # Aleft point is within B, collision.

    #print "unhandled!"
    return(False)

def sliding_window(listIn, window):
	# A cleaned up version of movingAverage()
	# Unlike movingAverage, only returns the y array this time.

    if window == 1: # just return the original array
        return(listIn)

    half_window_left = int(math.ceil(window / 2.0)) # correct for floating error division.
    half_window_right = int(math.floor(window / 2.0))

    y = []

    for n in xrange(half_window_left, len(listIn)-half_window_right):
        y.append(sum(listIn[n-half_window_left:n+half_window_right]))

    return(y)