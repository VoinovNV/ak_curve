cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(NAMESPACE ak_curve)
set(COMPONENT lib)
project(${NAMESPACE}_${COMPONENT}
    VERSION 0.0.0
    LANGUAGES C
)


add_library(${PROJECT_NAME}
    curve.h
    curve.c
)

add_library(${NAMESPACE}::${COMPONENT} ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")

