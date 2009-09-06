"""

gMaker.py

**Purpose**

Makes the actual wxWidgets GUI based on the gui defined by glbase

"""

# ----------------------------------------------------------------------
# list of gui definitions and types.
# ----------------------------------------------------------------------

valid_datatypes = ["int", "str", "loc"]

# ----------------------------------------------------------------------
# Entry point, call resulting_panel = make_gui(function)
# ----------------------------------------------------------------------
def make_gui(function, **kargs):
    """
    make the gui and return a formatted wxPanel ready for binding into the
    GUI

    **Arguments**

        function
            a function with a valid __gui__ attribute
    """
