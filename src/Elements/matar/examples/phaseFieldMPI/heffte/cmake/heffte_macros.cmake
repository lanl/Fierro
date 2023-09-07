# usage: heffte_find_libraries(REQUIRED foo1 foo2 OPTIONAL foo3 PREFIX bar LIST saloon NO_DEFAULT_PATH)
# this will search for foo1/2/3 in bar/<lib/lib64/arch> with NO_DEFAULT_PATH (skip if defaults are to be used)
# the found libraries will be added to a list heffte_saloon
# the search results will be also added to cached variables heffte_foo1 and heffte_foo2
# that can be used to call find_package_handle_standard_args(... DEFAULT_MSG heffte_foo1 heffte_foo2)
# the macro will not call find_package_handle_standard_args()
# the optional libraries will not create cached entries and mission optional will not be added to the LIST
macro(heffte_find_libraries)
    cmake_parse_arguments(heffte_findlibs "NO_DEFAULT_PATH" "PREFIX;LIST" "REQUIRED;OPTIONAL" ${ARGN})

    foreach(_tsg_lib ${heffte_findlibs_REQUIRED})
        list(APPEND heffte_findlibs_required "heffte_${_tsg_lib}")
    endforeach()

    set(heffte_findlibs_default "")
    if (${heffte_findlibs_NO_DEFAULT_PATH})
        set(heffte_findlibs_default "NO_DEFAULT_PATH")
    endif()

    foreach(_hfft_lib ${heffte_findlibs_REQUIRED} ${heffte_findlibs_OPTIONAL})

        find_library(heffte_${_hfft_lib} ${_hfft_lib}
                     HINTS "${heffte_findlibs_PREFIX}"
                     HINTS "${heffte_findlibs_PREFIX}/lib/"
                     HINTS "${heffte_findlibs_PREFIX}/lib/intel64"
                     HINTS "${heffte_findlibs_PREFIX}/${CMAKE_LIBRARY_ARCHITECTURE}/lib/"
                     HINTS "${heffte_findlibs_PREFIX}/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
                     HINTS "${heffte_findlibs_PREFIX}/lib64/"
                     HINTS "${heffte_findlibs_PREFIX}/${CMAKE_LIBRARY_ARCHITECTURE}/lib64/"
                     HINTS "${heffte_findlibs_PREFIX}/lib64/${CMAKE_LIBRARY_ARCHITECTURE}"
                     ${heffte_findlibs_default})

        if (CMAKE_FIND_DEBUG_MODE)
            message(STATUS "Heffte searching library: ${_hfft_lib} => ${heffte_${_hfft_lib}}")
        endif()

        list(FIND heffte_findlibs_required "heffte_${_hfft_lib}" heffte_findlibs_is_required)
        if (${heffte_findlibs_is_required} EQUAL -1) # not a required library
            if (heffte_${_hfft_lib})
                list(APPEND heffte_${heffte_findlibs_LIST} ${heffte_${_hfft_lib}})
            endif()
            unset(heffte_${_hfft_lib} CACHE)
        else()
            list(APPEND heffte_${heffte_findlibs_LIST} ${heffte_${_hfft_lib}})
        endif()

    endforeach()

    foreach(_hfft_lib default required NO_DEFAULT_PATH PREFIX LIST REQUIRED OPTIONAL) # cleanup
        unset(heffte_${_hfft_lib})
    endforeach()
    unset(_hfft_lib)
endmacro()
