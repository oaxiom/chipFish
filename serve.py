"""
A web interface for chipfish.

This is the abandonment of the GUI

"""
import time, os, sys
import cli as chipFish_draw
from glbase_wrapper import location
from wsgiref.simple_server import make_server

HOST_NAME = ""
PORT_NUMBER = 8081

cf = chipFish_draw.cfApp()

pythonpath = os.getenv("PYTHONPATH", '')
pythonpath += ":" + os.path.abspath(sys.path[0]).split()[0]
os.environ["PYTHONPATH"] = pythonpath
os.chdir("htdocs/") # The location I want to serve from

class yApp:
    """Produce the same output, but using a class
    """
    def __init__(self, environ, start_response):
        self.environ = environ
        self.start = start_response

    def __iter__(self):
        self.environ["PATH"] = "/Users/hutchinsandrew/chipfish/htdocs/"
        self.environ["PWD"] = "/Users/hutchinsandrew/chipfish/htdocs/"
        self.environ["HOME"] = "/Users/hutchinsandrew/chipfish/htdocs/"
        os.chdir("/Users/hutchinsandrew/chipfish/htdocs/")
        print self.environ
        status = '200 OK'
        response_headers = [('Content-type','text/html')]
        self.start(status, response_headers)
        cf.draw.setLocation(loc=location(loc="chr19:34010332-34653006"))
        #cf.draw.exportImage("/Users/hutchinsandrew/chipfish/htdocs/last_image.png")
        ret = """
        <html><head><title>chipFish</title></head>\n
        <body><p>This is a test...</p>\n
        <p><img src='/Users/hutchinsandrew/chipfish/htdocs/tmp/last_image.png' alt='browser image genome' /></p>\n

        </body></html>\n
        """
        yield(ret)

httpd = make_server('', 8080, yApp)
print "Serving HTTP on port 8080..."

# Respond to requests until process is killed
httpd.serve_forever()


"""



import BaseHTTPServer
import SimpleHTTPServer

import cli as chipFish_draw
from glbase_wrapper import location

HOST_NAME = ""
PORT_NUMBER = 8081

cf = chipFish_draw.cfApp()

pythonpath = os.getenv("PYTHONPATH", '')
pythonpath += ":" + os.path.abspath(sys.path[0]).split()[0]
os.environ["PYTHONPATH"] = pythonpath
os.chdir("htdocs/") # The location I want to serve from

class MyHandler(BaseHTTPServer.BaseHTTPRequestHandler):
    def do_HEAD(s):
        s.send_response(200)
        s.send_header("Content-type", "text/html")
        s.end_headers()

    def do_GET(self):
        #Respond to a GET request.
        #self.path = os.path.join(os.path.expanduser("~"), "chipFish/htdocs/")

        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        self.wfile.write("<html><head><title>chipFish</title></head>\n")
        self.wfile.write("<body><p>This is a test...</p>\n")
        self.wfile.write("<p>You accessed path: %s</p>\n" % self.path)

        #if self.a and not self.a:

        cf.draw.setLocation(loc=location(loc="chr19:34010332-34653006"))
        #cf.draw.exportImage("tmp/last_image.png")
        self.wfile.write("<p><img src='tmp/last_image.png' alt='browser image genome' /></p>\n")

        self.wfile.write("</body></html>\n")


if __name__ == '__main__':
    server_class = BaseHTTPServer.HTTPServer
    httpd = server_class((HOST_NAME, PORT_NUMBER), MyHandler)
    print time.asctime(), "Server Starts - %s:%s" % (HOST_NAME, PORT_NUMBER)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    httpd.server_close()
    print time.asctime(), "Server Stops - %s:%s" % (HOST_NAME, PORT_NUMBER)



def httpd_serve_forever(port=8080) :
    #Start a webserver on a local port.
    import BaseHTTPServer
    import CGIHTTPServer

    class __HTTPRequestHandler(CGIHTTPServer.CGIHTTPRequestHandler):
        def is_cgi(self) :
            if self.path == "/create.cgi":
                self.cgi_info = '', 'create.cgi'
                return True
            return False

    # Add current directory to PYTHONPATH. This is
    # so that we can run the standalone server
    # without having to run the install script.
    pythonpath = os.getenv("PYTHONPATH", '')
    pythonpath += ":" + os.path.abspath(sys.path[0]).split()[0]
    os.environ["PYTHONPATH"] = pythonpath

    htdocs = "htdocs/path/to/"
    os.chdir(htdocs)

    HandlerClass = __HTTPRequestHandler
    ServerClass = BaseHTTPServer.HTTPServer
    httpd = ServerClass(('', port), HandlerClass)
    print "Serving HTTP on localhost:%d ..." % port

    try :
        httpd.serve_forever()
    except KeyboardInterrupt:
        sys.exit(0)
# end httpd_serve_forever()
""" # From weblog.
