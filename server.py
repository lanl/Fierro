import os
import cherrypy

PATH = os.path.abspath(os.path.dirname(__file__))
class Root(object): pass

cherrypy.tree.mount(Root(), '/', config={
    '/': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': PATH,
            'tools.staticdir.index': 'index.html',
        },
})

cherrypy.engine.start()