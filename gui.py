"""
gui, part of chipFish

(c) 2008-2009 oAxiom

Not for distribution.

Class container for constants and other data

"""

import wx, guiWxParent, math, sys

import ctypes
import cairo
from ctypes.util import find_library

cairo_dll = ctypes.CDLL(find_library("cairo"))

# Pycairo's API representation (from pycairo.h)
class Pycairo_CAPI(ctypes.Structure):
   _fields_ = [
      ('Context_Type', ctypes.py_object),
      ('Context_FromContext', ctypes.PYFUNCTYPE(ctypes.py_object,
                                                ctypes.c_void_p,
                                                ctypes.py_object,
                                                ctypes.py_object)),
      ('FontFace_Type', ctypes.py_object),
      ('FontFace_FromFontFace', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('FontOptions_Type', ctypes.py_object),
      ('FontOptions_FromFontOptions', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Matrix_Type', ctypes.py_object),
      ('Matrix_FromMatrix', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Path_Type', ctypes.py_object),
      ('Path_FromPath', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Pattern_Type', ctypes.py_object),
      ('SolidPattern_Type', ctypes.py_object),
      ('SurfacePattern_Type', ctypes.py_object),
      ('Gradient_Type', ctypes.py_object),
      ('LinearGradient_Type', ctypes.py_object),
      ('RadialGradient_Type', ctypes.py_object),
      ('Pattern_FromPattern', ctypes.c_void_p),
      ('ScaledFont_Type', ctypes.py_object),
      ('ScaledFont_FromScaledFont', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Surface_Type', ctypes.py_object),
      ('ImageSurface_Type', ctypes.py_object),
      ('PDFSurface_Type', ctypes.py_object),
      ('PSSurface_Type', ctypes.py_object),
      ('SVGSurface_Type', ctypes.py_object),
      ('Win32Surface_Type', ctypes.py_object),
      ('XlibSurface_Type', ctypes.py_object),
      ('Surface_FromSurface', ctypes.PYFUNCTYPE(ctypes.py_object, ctypes.c_void_p)),
      ('Check_Status', ctypes.PYFUNCTYPE(ctypes.c_int, ctypes.c_int))]

# look up the API
ctypes.pythonapi.PyCObject_Import.restype = ctypes.POINTER(Pycairo_CAPI)
pycairo_api = ctypes.pythonapi.PyCObject_Import("cairo", "CAPI").contents;

ContextType = pycairo_api.Context_Type

def Context_FromSWIGObject(swigObj):
    ptr = ctypes.c_void_p(int(swigObj))
    #increment the native context's ref count, since the Pycairo_Context decrements it
    #when it is finalised.
    cairo_dll.cairo_reference(ptr)
    return pycairo_api.Context_FromContext(ptr, ContextType, None)

class drawPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, size=(-1,-1), style=wx.FULL_REPAINT_ON_RESIZE) # This has to be set to full repaint :)
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        # put it into a sizer to make it movable.
        #sizer = wx.BoxSizer(wx.VERTICAL)
        #sizer.Add(self.canvas, 1, wx.LEFT|wx.TOP|wx.GROW)
        #self.SetSizer(sizer)
        #self.SetSizer(sizer)
        self.Fit()

    def OnPaint(self, event):
        # override;
        dc = wx.PaintDC(self) # make each time OnPaint is called.
        #dc = wx.AutoBufferedPaintDC(self) # not clear why this doesn't work ...
        w,h = dc.GetSizeTuple()
        gc = wx.GraphicsContext.Create(dc)
        nc = gc.GetNativeContext()
        ctx = Context_FromSWIGObject(nc)

        # Cairo:
        ctx.set_source_rgba(0,0,0,1)
        ctx.rectangle(0,0,w,h)
        ctx.set_source_rgba(1, 1, 0, 0.80)
        ctx.rectangle(0, 0, 0.5, 0.5)
        ctx.fill()
        ctx.arc(2.*w/3,2.*h/3.,min(w,h)/4. - 10,0, math.pi*2)
        ctx.set_source_rgba(0,1,1,0.5)
        ctx.fill_preserve()
        ctx.set_source_rgb(1,0.5,0)
        #ctx.stroke()

        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(0, 0)
        ctx.line_to(w, h)
        ctx.move_to(1, 0)
        ctx.line_to(0, 1)
        ctx.set_line_width(0.2)
        ctx.stroke()

        ctx.rectangle(0, 0, 0.5, 0.5)
        ctx.set_source_rgba(1, 0, 0, 0.80)
        ctx.fill()

        ctx.rectangle(0, 0.5, 0.5, 0.5)
        ctx.set_source_rgba(0, 1, 0, 0.60)
        ctx.fill()

        ctx.rectangle(0.5, 0, 0.5, 0.5)
        ctx.set_source_rgba(0, 0, 1, 0.40)
        ctx.fill()
        #self.Layout()

    def OnRedraw(self,evt):
        self.OnPaint(evt)

    def onEraseBackground(self, evt):
        # this is supposed to prevent redraw flicker on some X servers...
        # Still required?
        pass

class mainFrame(guiWxParent.frame_mainFrame_parent):
    def __init__(self, *args, **kwds):
        """
        The mainFrame for the app.
        """
        #sys.stderr = self # silence errors.
        #sys.stdout = self
        guiWxParent.frame_mainFrame_parent.__init__(self, None, -1) # inherit
        self.Maximize(True)
        # set-up the gDraw Cairo panel
        sizer = wx.BoxSizer(wx.VERTICAL)
        self._Drawer = drawPanel(self.gDrawPanel)
        sizer.Add(self._Drawer, 2, wx.EXPAND, 0)
        self.gDrawPanel.SetSizer(sizer)
        self.Fit()
        #
        #gDrawPanelSizer.Add(self.gDrawPanel, 2, wx.EXPAND, 0)

if __name__ == "__main__":
    app = wx.App()
    f = mainFrame(None, -1, "")
    f.Show()
    app.MainLoop()
