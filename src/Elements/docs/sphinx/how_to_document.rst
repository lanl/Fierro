How to add to docs
==================

How are the docs generated?
---------------------------
The process of generating the ELEMENTS documentation is based on `this blog <https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/>`_.

Why do it this way?
-------------------
The reasoning for why we produce our documentation in such a complicated way goes as follows. 
Sphinx is our preferred tool for documentation: 

1. over Doxygen (by itself) because Doxygen conduces developers to write their documentation in the form of an API reference manual, often leading to documentation that exhaustively describes the components of the software but, importantly, not how a user can compose them in order to do what he wants;
2. over Markdown-based documentation tools because Markdown is missing features that reStructuredText (in Sphinx) has.

Our hope is that with Sphinx, we can use the extended feature set of reStructuredText to write effective documentation that not only describes the parts of our software but also explains how to put them together. 
The complication is that Sphinx cannot parse C++ source directly, so we must use Doxygen to parse the C++ source into XML files and Breathe, an extension of Sphinx, to parse those XML files so that they can be used by Sphinx.

How to add a page like this one?
--------------------------------
1. Create a new reStructuredText file in the ``docs/sphinx/source`` directory.
2. Add the file name to the list of files under ``Contents:`` in the ``docs/sphinx/source/index.rst`` file.
3. Add the file to the list of files tracked in the ``docs/CMakeLists.txt`` file (beneath the ``DEPENDS`` keyword in the ``add_custom_command`` function in the Sphinx section of the file).
4. Write the page using reStructuredText.

How to write a page of reStructuredText?
-----------------------------------------
Please see the reStructuredText `primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_ on the Sphinx web site.
Also, when editing the reStructuredText source for this project's documentation, please:

* put each sentence on a new line, to keep the "diffs" sensible; and
* use a spell checker.

How to include source documentation?
------------------------------------
1. Write source documentation in the headers using comments and comments blocks formatted for Doxygen, as explained on `this page <https://www.doxygen.nl/manual/docblocks.html>`_ in the Doxygen documentation.
2. Add `Breathe's Doxygen directives <https://breathe.readthedocs.io/en/latest/directives.html>`_, e.g. for documented namespaces, classes, structs, functions, to the Sphinx reStructuredText files.

How to generate documentation locally?
--------------------------------------
1. Ensure that documentation generation dependencies are installed.
   Recall that these are `Doxygen <https://doxygen.nl/download.html>`_, `Sphinx <https://www.sphinx-doc.org/en/master/usage/installation.html>`_, `Breathe <https://breathe.readthedocs.io/en/latest/#download>`_, and `Read the Docs Sphinx Theme <https://sphinx-rtd-theme.readthedocs.io/en/stable/installing.html>`_. 
2. Configure your CMake build of ELEMENTS with the option flags ``-DWITH_DOCS=ON`` and, if you plan to install the documentation somewhere, ``DCMAKE_INSTALL_DOCDIR=/path/to/your/docs/install``.
3. Build ELEMENTS.

How to publish changes to GitHub Pages?
---------------------------------------
This is done automatically.
When commits are pushed to the ``master`` branch, several GitHub Actions workflows are launched.
One of these is a workflow that regenerates the documentation, commits, and pushes the updates to the ``gh-pages`` branch, which is sourced by GitHub Pages.

