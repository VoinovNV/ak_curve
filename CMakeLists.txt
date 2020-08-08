cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(ak_curve
    VERSION 0.0.0
    LANGUAGES C
)

set(SUBPROJECTS
    ak_curve_lib
    ak_curve_tests
    ak_curve_tests_check
    ak_curve_demo
)

macro(find_package)
    if(NOT "${ARGV0}" IN_LIST SUBPROJECTS)
        _find_package(${ARGV})
    endif()
endmacro()

foreach(subdir IN LISTS SUBPROJECTS)
    add_subdirectory(${subdir})
endforeach()
