find_package(Doxygen REQUIRED OPTIONAL_COMPONENTS dot)

set(DOXYGEN_GENERATE_HTML      "YES")
set(DOXYGEN_DISABLE_INDEX       "NO")
set(DOXYGEN_GENERATE_TREEVIEW  "YES") # left-hand index
set(DOXYGEN_SORT_GROUP_NAMES   "YES")
set(DOXYGEN_SORT_BRIEF_DOCS     "NO")
set(DOXYGEN_SORT_MEMBER_DOCS    "NO")
set(DOXYGEN_COLLABORATION_GRAPH "NO")
set(DOXYGEN_INCLUDE_GRAPH       "NO")
set(DOXYGEN_INCLUDED_BY_GRAPH   "NO")
set(DOXYGEN_GRAPHICAL_HIERARCHY "NO")
set(DOXYGEN_CLASS_GRAPH         "NO")
set(DOXYGEN_CLASS_DIAGRAMS      "NO")

set(DOXYGEN_HTML_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/../html/")

set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")
set(DOXYGEN_HTML_EXTRA_STYLESHEET  "${CMAKE_CURRENT_SOURCE_DIR}/heffte.css")

set(DOXYGEN_HTML_COLORSTYLE_HUE    "205")
set(DOXYGEN_HTML_COLORSTYLE_SAT     "95")
set(DOXYGEN_HTML_COLORSTYLE_GAMMA  "190")

# PREDEFINED is used only during doxygen build to show all methods and classes for all options
set(DOXYGEN_PREDEFINED "Heffte_ENABLE_FFTW;Heffte_ENABLE_MKL;Heffte_ENABLE_CUDA;Heffte_ENABLE_ROCM;Heffte_ENABLE_ONEAPI;Heffte_ENABLE_TRACING")

# better HeFFTe name
set(DOXYGEN_PROJECT_NAME "Highly Efficient FFT for Exascale: HeFFTe v${Heffte_VERSION_MAJOR}.${Heffte_VERSION_MINOR}")
set(DOXYGEN_PROJECT_NUMBER "")

doxygen_add_docs(heffte_doxygen
                 README.md
                 doxygen/installation.md
                 doxygen/basic_usage.md
                 python/README.md
                 fortran/README.md
                 include/
                 src/heffte_plan_logic.cpp
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../
                 COMMENT "Building the ${PROJECT_NAME} documentation")

add_dependencies(Heffte heffte_doxygen)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../html DESTINATION "share/heffte/docs")
