cmake_minimum_required(VERSION 3.16.3)
project(curvature_optimizer)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()

set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

get_directory_property(PARENT_DIR PARENT_DIRECTORY)
if(NOT PARENT_DIR)
    set(TOPLEVEL_PROJECT ON)
else()
    set(TOPLEVEL_PROJECT OFF)
endif()

option(BUILD_HEADLESS "Build without GUI" OFF)
option(BUILD_TESTS "Build with tests" ON)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

if (EMSCRIPTEN)
    set(BUILD_HEADLESS OFF CACHE BOOL "Disabled in Emscripten build" FORCE)
    set(BUILD_TESTS OFF CACHE BOOL "Disabled in Emscripten build" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Disabled in Emscripten build" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    if(WIN32)
        if (MSVC)
            add_compile_options(/bigobj)
        else ()
            add_compile_options(-Wa,-mbig-obj)
        endif ()
    endif()
    if((UNIX OR APPLE) AND CMAKE_COMPILER_IS_GNUCXX)
        add_compile_options(-ggdb)
    endif()
endif()

add_subdirectory(deps)
add_subdirectory(src)

include(clang-format)