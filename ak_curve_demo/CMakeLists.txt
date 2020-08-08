cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(NAMESPACE ak_curve)
set(COMPONENT demo)
project(${NAMESPACE}_${COMPONENT}
    VERSION 0.0.0
    LANGUAGES C
)

add_executable(${PROJECT_NAME}
    Source.c
)

add_executable(${NAMESPACE}::${COMPONENT} ALIAS ${PROJECT_NAME})

find_package(ak_curve_lib 0.0.0 REQUIRED)
find_package(benchmark REQUIRED)
target_link_libraries(${PROJECT_NAME} ak_curve::lib benchmark::benchmark_main)

find_library(LIBAKR libakrypt-static.a
    PATHS /home/user/Projects/ak_curve/libakrypt-0.x/
)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../libakrypt-0.x/source/>)
target_link_libraries(${PROJECT_NAME} ${LIBAKR})

