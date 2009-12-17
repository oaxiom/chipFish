"""
some choice utils from glbase

"""

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
