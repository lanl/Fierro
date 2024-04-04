#!/bin/bash -e

# Check if the uncrustify build directory exists and is not empty in the parent directory; if not, clone it

if [ ! -d "${UNCRUSTIFY_SOURCE_DIR}" ];
then
    echo "Missing Uncrustify submodule, downloading...."
    git submodule update ${UNCRUSTIFY_SOURCE_DIR}
fi

if [ ! -f "${UNCRUSTIFY_BUILD_DIR}/uncrustify" ]; then
    echo "Uncrustify not built, installing..."
    cmake -D CMAKE_BUILD_TYPE=Release -B "${UNCRUSTIFY_BUILD_DIR}" -S "${UNCRUSTIFY_SOURCE_DIR}"
    make -C "${UNCRUSTIFY_BUILD_DIR}" -j${FIERRO_BUILD_CORES}
else
    echo "Uncrustify already exists and has been built"
fi
