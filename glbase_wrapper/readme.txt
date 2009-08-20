
This is a part of chipFish, it is a series of headers overriding the vanilla glbase.

These are not released under the GPL, and as they inherit from glbase, they do not need to be. 

(c) 2009 A.Hutchins and oAxiom.

Introduction
------------

This is a series of wrappers around the vanilla glbase that add on things like gui definitions and localised documentation strings. These are then plugged into chipFish. which knows how to interpret them.

Eventually here I'll put a definition of the gui (it remains in flux at the current writing).

GUI Specifications:
-------------------

object.method.__doc__ = "string"

Override the glbase documentation and insert a reformatted explanation that goes into the interface of glbase

object.method.__tooltip__ = "string"

the tooltip describing in a single pithy sentence the rols of this function.



