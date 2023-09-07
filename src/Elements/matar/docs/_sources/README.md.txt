## Requirements

1. doxygen -- Built with version 1.9.4 
2. sphinx -- Built with version 4.5.0

## Running
To run the documentation program in the docs_sphinx run `sh docify.sh` this does the following 

    1. Runs doxygen which generates html and xml files used for the website 
    
    2. Runs sphinx which primarily adds style headers and different page naviagation to the site 
    
    3. Moves the output from `docs_sphinx/_build/html` to `docs/` which is where github is looking for a webiste 
    
    4. Moves the output from `docs_sphinx/_build/html/_static` to `docs/static`. For some reason github can't find the folder `_static` but can find `static`
    
    5. Renames the location of .css files on each page from `_static/` to static/ to account for this. 

## Known issues

    The second step, running sphinx, generates a lot of errors. 
 
    Many of the default links in particular "view source code" is broken. But all the links to classes/macros/functions should work

    The above procudure in the previous section is really awkward and there is probably a way in sphinx to specify the output directory and what to name the `_static` libary.
