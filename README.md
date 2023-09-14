# Pages

This branch holds the HTML, CSS, and JS that controls the webpage. Modifications to this branch will automatically be reflected in the webpage after a short delay for redeployment. 

## Editing this content

When making changes to this page, you will most likely find yourself editing the text in one of the HTML files, editing the CSS in assets/css/style.css, and/or adding images under assets/images/. If there is already an example of what you want on the page, you are encouraged to copy and modify it for your own content.

### Local Development
To facilitate local development, there is a python file `server.py` that will host the web page locally so you can see the changes you make in real-time as well as use your browsers built-in development tools. Most changes you make should be tested/inspected by using this server before pushing to the remote branch.

To run the development server, first ensure you have some version of python installed (typically provided through Anaconda on Windows). Then install the `cherrypy` dependency.

```pip install cherrypy```
or 
```conda install cherrypy```

Then you can run the server with 

```python server.py```

To access the locally hosted web page, open a browser and navigate to the server host location (typically `http://127.0.0.1:8080/`). To ensure that your changes are consistently reflected on browser refresh, you should disable caching under the network tab of the browsers dev tools window accessed through f12.