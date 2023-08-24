# The process of generating documentation here is based on a blog post titled
# "Clear, Functional C++ Documentation with Sphinx + Breathe + Doxygen + CMake"
# (https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/).

# The reasoning for why we produce our documentation in such a complicated way
# goes as follows. Sphinx is our preferred tool for documentation: 
# (1) over Doxygen (by itself) because Doxygen conduces developers to write
# their documentation in the form of an API reference manual, often leading to
# documentation that exhaustively describes the components of the software but,
# importantly, not how a user can compose them in order to do what he wants;
# (2) over Markdown-based documentation tools because Markdown is missing
# features that reStructuredText (in Sphinx) has.
# Our hope is that with Sphinx, we can use the extended feature set of
# reStructuredText to write effective documentation that not only describes the
# parts of our software but also explains how to put them together. The
# complication is that Sphinx cannot parse C++ source directly, so we must use
# Doxygen to parse the C++ source into XML files and Breathe, an extension of
# Sphinx, to parse those XML files so that they can be used by Sphinx.


################################################################################
# Doxygen
################################################################################
# List the header files to be included in the documentation
file(GLOB_RECURSE DOC_HEADERS ../include/*.h)

# Set up variables for Doxygen template configuration file, configuration file,
# and output html file
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile) 
set(DOXYGEN_INPUT_DIR "${PROJECT_SOURCE_DIR}/common ${PROJECT_SOURCE_DIR}/elements ${PROJECT_SOURCE_DIR}/geometry ${PROJECT_SOURCE_DIR}/slam ${PROJECT_SOURCE_DIR}/swage")
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)

# Replace the placeholders in the template configuration file template with the
# values of the variables defined above
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Make the output directory (Doxygen won't do so automatically)
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}/doxygen)

# Add rule for running the Doxygen command. Only regenerate Doxygen when the
# Doxyfile or documented headers change
add_custom_command(
  OUTPUT ${DOXYGEN_INDEX_FILE}
  DEPENDS ${DOC_HEADERS}
  COMMAND Doxygen::doxygen ${DOXYFILE_OUT}
  MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
  COMMENT "Processing header documentation with Doxygen"
  VERBATIM
)

# Add target that ensures the rule above is applied (command is executed)
add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})


################################################################################
# Sphinx with Breathe
################################################################################
# Set up variables for Sphinx documentation generation
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

# Add rule for running the Sphinx command. Only regenerate Sphinx when:
# - Doxygen has rerun;
# - the doc files have been updated; or
# - the Sphinx config has been updated.
add_custom_command(
  OUTPUT ${SPHINX_INDEX_FILE}
  COMMAND
    ${SPHINX_EXECUTABLE} -b html
    # Tell Breathe where to find the Doxygen output
    -Dbreathe_projects.ELEMENTS=${DOXYGEN_OUTPUT_DIR}/xml
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    # Other docs files you want to track should go here (or in some variable)
    ${SPHINX_SOURCE}/index.rst
    ${SPHINX_SOURCE}/how_to_install.rst
    ${SPHINX_SOURCE}/how_to_use.rst
    ${SPHINX_SOURCE}/elements.rst
    ${SPHINX_SOURCE}/geometry.rst
    ${SPHINX_SOURCE}/matar.rst
    ${SPHINX_SOURCE}/slam.rst
    ${SPHINX_SOURCE}/swage.rst
    ${SPHINX_SOURCE}/examples.rst
    ${SPHINX_SOURCE}/how_to_document.rst
    ${DOXYGEN_INDEX_FILE}
  MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
  COMMENT "Generating documentation with Sphinx"
)

# Add target that ensures the rule above is applied (command is executed)
add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

# Install the Sphinx documentation where specified by configuration
install(DIRECTORY ${SPHINX_BUILD}/ DESTINATION ${CMAKE_INSTALL_DOCDIR})
