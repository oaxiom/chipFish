"""
A web interface for chipfish.

Copyright 2011-2014 Andrew Hutchins

This is the abandonment of the GUI

BUGS:
-----
. If the local machine is not connected to the internet, it sometimes takes an age to begin serving
. favicon.ico not found
. root dir is not htdocs, chipfish .
. This must surely be a really bad way to do this.
. Isn't there a better (more secure) lightweight way to do this?
. GUI? Even possible? - not with this crappy server!

"""
import time, os, sys
from app import app
from glbase_wrapper import location

if len(sys.argv) < 1:
    print("serve.py <track_list_file.txt>")
    quit()

cf = app() # Startup chipFish
cf.startup(sys.argv[1])

__version__ = "0.1"

__all__ = ["HTTPRequestHandler"]

import os, sys, time
import posixpath
import http.server
import urllib.request, urllib.parse, urllib.error
import cgi
import shutil
import mimetypes

try:
    from io import StringIO
except ImportError:
    from io import StringIO

class chipfishHTTPServerHandler(http.server.BaseHTTPRequestHandler):
    """
    HTTP request handler with GET and HEAD commands.

    The GET and HEAD requests are identical except that the HEAD
    request omits the actual contents of the file.
    """
    server_version = "chipfishHTTPServer/" + __version__

    def do_GET(self):
        """
        Serve a GET request.
        """
        f = self.send_head()
        if f:
            self.copyfile(f, self.wfile)
            f.close()

    def do_HEAD(self):
        """
        Serve a HEAD request.
        """
        f = self.send_head()
        if f:
            f.close()

    def send_head(self):
        """
        Common code for GET and HEAD commands.

        This sends the response code and MIME headers.

        Return value is either a file object (which has to be copied
        to the outputfile by the caller unless the command was HEAD,
        and must be closed by the caller under all circumstances), or
        None, in which case the caller has nothing further to do.

        """
        path, query = self.translate_get(self.path)

        print(path, query)

        if query:
            if "loc" in query and query["loc"]:
                self.update_app(query["loc"])
            elif "chr" in query and "left" in query and "right" in query:
                self.update_app(loc=location(chr=query["chr"], left=query["left"], right=query["right"]))

        f = None
        if os.path.isdir(path):
            if not self.path.endswith('/'):
                # redirect browser - doing basically what apache does
                self.send_response(301)
                self.send_header("Location", self.path + "/")
                self.end_headers()
                return(None)

            for index in "index.html", "index.htm":
                index = os.path.join(path, index)
                if os.path.exists(index):
                    path = index
                    break
            else: # No index, serve a directory listing
                self.send_error(404, "File not found")

        ctype = self.guess_type(path)
        try:
            # Always read in binary mode. Opening files in text mode may cause
            # newline translations, making the actual size of the content
            # transmitted *less* than the content-length!
            f = open(path, 'rb')
        except IOError:
            self.send_error(404, "File not found")
            return None

        self.send_response(200)
        self.send_header("Content-type", ctype)
        fs = os.fstat(f.fileno())
        self.send_header("Content-Length", str(fs[6]))
        self.send_header("Last-Modified", self.date_time_string(fs.st_mtime))
        self.end_headers()
        return f

    def update_app(self, loc=None):
        assert loc, "no location!"

        cf.draw.setLocation(loc=location(loc=loc))
        cf.draw.exportImage(os.path.join(sys.path[0], "htdocs/tmp/last_image.png"))

    def translate_get(self, get):
        """Translate a /-separated PATH to the local filename syntax.

        Components that mean special things to the local file system
        (e.g. drive or directory names) are ignored.  (XXX They should
        probably be diagnosed.)

        """
        path = get.split('?',1)[0]
        query = None
        if "?" in get:
            cmds = get.split("?")[1].split("&")
            query = {q.split("=")[0]: q.split("=")[1] for q in cmds}
        path = path.split('#',1)[0]
        path = posixpath.normpath(urllib.parse.unquote(path))
        words = path.split('/')
        words = [_f for _f in words if _f]
        path = os.getcwd()
        for word in words:
            drive, word = os.path.splitdrive(word)
            head, word = os.path.split(word)
            if word in (os.curdir, os.pardir): continue
            path = os.path.join(path, word)

        return (path, query)

    def copyfile(self, source, outputfile):
        """Copy all data between two file objects.

        The SOURCE argument is a file object open for reading
        (or anything with a read() method) and the DESTINATION
        argument is a file object open for writing (or
        anything with a write() method).

        The only reason for overriding this would be to change
        the block size or perhaps to replace newlines by CRLF
        -- note however that this the default server uses this
        to copy binary data as well.

        """
        shutil.copyfileobj(source, outputfile)

    def guess_type(self, path):
        """Guess the type of a file.

        Argument is a PATH (a filename).

        Return value is a string of the form type/subtype,
        usable for a MIME Content-type header.

        The default implementation looks the file's extension
        up in the table self.extensions_map, using application/octet-stream
        as a default; however it would be permissible (if
        slow) to look inside the data to make a better guess.

        """

        base, ext = posixpath.splitext(path)
        if ext in self.extensions_map:
            return self.extensions_map[ext]
        ext = ext.lower()
        if ext in self.extensions_map:
            return self.extensions_map[ext]
        else:
            return self.extensions_map['']

    if not mimetypes.inited:
        mimetypes.init() # try to read system mime.types
    extensions_map = mimetypes.types_map.copy()
    extensions_map.update({
        '': 'application/octet-stream', # Default
        '.py': 'text/plain',
        '.c': 'text/plain',
        '.h': 'text/plain',
        })

def serve(handler_class=chipfishHTTPServerHandler, server_class=http.server.HTTPServer):
    # Hack the port.
    sys.argv[1] = 8080
    
    server_address = ('127.0.0.1', 8080)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()

if __name__ == '__main__':
    serve()
