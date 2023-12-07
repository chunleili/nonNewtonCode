cmake_minimum_required(VERSION 3.0)

project(GenericParameters)

set(PROJECT_PATH ${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(${PROJECT_SOURCE_DIR}/CMake/Common.cmake)

add_subdirectory(GenericParameters)
