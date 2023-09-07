How to install
==============
These instructions are written for installation from the terminal of a Unix-like (MacOS, Linux, etc.) system.


Clone
-----
To clone the **ELEMENTS** repository, enter the following command. ::

  git clone --recursive https://github.com/lanl/ELEMENTS.git

The ``--recursive`` flag ensures that the submodules are initialized as the repository is cloned.
If you forget to add this flag, you can initialize the submodules later by entering the following command. ::

  git submodule update --init

Then enter the ELEMENTS directory. ::

  cd ELEMENTS


Configure
---------
It is often convenient to write a configuration script to avoid having to write out option flags more than once and to have a record of what options were used in the build.
A starting point for a configuration script might be::

  #!/bin/bash
  ELEMENTS_DIR=/path/to/elements/repo
  cmake \
    ${ELEMENTS_DIR}

Suppose you write this to a file called ``my_config.sh``.
To make the script executable, run the following command. ::

  chmod +x my_config.sh

Configure to build in place
^^^^^^^^^^^^^^^^^^^^^^^^^^^
To build **ELEMENTS** in place, substitute the placeholder ``/path/to/elements/repo`` in the configuration script for ``.``, specifying the current working directory. 
Then simply run the script as follows. ::

  ./my_config.sh

Configure to build elsewhere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To build **ELEMENTS**, say, in a ``build`` subdirectory, enter the following commands. ::
  
  mkdir build; cd build

Substitute the placeholder ``/path/to/elements/repo`` in the configuration script for ``..``, specifying the directory above the current working directory. 
Copy the configuration script into the build directory, enter it, and run the configuration. ::

  cp my_config.sh build; cd build; ./my_config.sh

Configure to install in place
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Whether building in place or in a separate build directory, you add configure CMake to install in place.
The option flag for this is ``-DCMAKE_INSTALL_PREFIX=\`pwd\```.
Once added to the configuration script, it becomes the following. ::

  #!/bin/bash
  ELEMENTS_DIR=/path/to/elements/repo
  cmake \
    -DCMAKE_INSTALL_PREFIX=`pwd` \
    ${ELEMENTS_DIR}

Enable optional dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable the VTK dependency, used for reading and writing high-order hexahedral meshes in the VTK format, add the following two options to the configuration script: ``-DWITH_VTK=ON`` and ``-DVTK_DIR=/path/to/vtk/install``.
In the second option, you should substitute the placeholder ``/path/to/vtk/install`` with the actual path to your VTK installation.

To enable the use of BLAS/LAPACK in ELEMENTS, add the following option to the configuration script: ``-DWITH_BLAS_LAPACK=ON``.

To build the documentation locally, add the following option to the configuration script: ``-DWITH_DOCS=ON``. 
Bear in mind that there are 3 dependencies involved in generating the documentation: Doxygen, Sphinx, and Breathe.
The easiest way to install Doxygen is using whatever package manager you have on your OS, if it includes Doxygen.
Since Sphinx and Breathe are Python-based, the easiest way to install them is using ``pip`` specifically.


Build & Install
---------------
To build **ELEMENTS**, first configure as explained above and then enter ``make`` in the build directory.
To install the **ELEMENTS** libraries, binaries, and header files in the directory specified by the configuration, enter ``make install`` in the build directory.
