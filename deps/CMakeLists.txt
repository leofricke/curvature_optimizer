include(FetchContent)

if(EMSCRIPTEN)
    add_compile_options(-pthread)
endif()

if(NOT TARGET eigen)
    message(STATUS "Adding eigen dependency")

    FetchContent_Declare(
        eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz
        SYSTEM
        EXCLUDE_FROM_ALL
        SOURCE_SUBDIR thisdoesnotexist
    )
    FetchContent_MakeAvailable(eigen)

    add_library(eigen INTERFACE EXCLUDE_FROM_ALL)
    target_include_directories(eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS eigen)
    
    set(EIGEN3_INCLUDE_DIR ${eigen_SOURCE_DIR}) # make pmp use this version of eigen
endif()

if(NOT TARGET pmp)
    message(STATUS "Adding pmp-library dependency")

    set(PMP_BUILD_VIEWERS NOT BUILD_HEADLESS CACHE BOOL "" FORCE)
    set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(PMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(PMP_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(PMP_INSTALL OFF CACHE BOOL "" FORCE)
    set(PMP_STRICT_COMPILATION OFF CACHE BOOL "" FORCE)

    add_compile_definitions(TUDO)

    FetchContent_Declare(
        pmp-library
        GIT_REPOSITORY https://github.com/leofricke/pmp-library
        GIT_TAG main
        SYSTEM
        EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(pmp-library)
endif()

if(NOT BUILD_HEADLESS AND NOT TARGET implot)
    message(STATUS "Adding implot dependency")

    add_subdirectory(implot SYSTEM EXCLUDE_FROM_ALL)
endif()

if(NOT TARGET nlohmann_json)
    message(STATUS "Adding json dependency")

    add_subdirectory(json SYSTEM EXCLUDE_FROM_ALL)
endif()

if(NOT TARGET happly)
    message(STATUS "Adding happly dependency")

    add_subdirectory(happly SYSTEM EXCLUDE_FROM_ALL)
endif()

if(NOT TARGET args)
    message(STATUS "Adding args dependency")

    add_subdirectory(args SYSTEM EXCLUDE_FROM_ALL)
endif()

if(BUILD_TESTS AND NOT TARGET GTest::gtest_main)
    message(STATUS "Adding googletest dependency")

    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/releases/download/v1.18.0/googletest-1.18.0.tar.gz
        SYSTEM
        EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(googletest)

    set_target_properties(gtest
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY     "${CMAKE_BINARY_DIR}/bin"
        LIBRARY_OUTPUT_DIRECTORY     "${CMAKE_BINARY_DIR}/lib"
        ARCHIVE_OUTPUT_DIRECTORY     "${CMAKE_BINARY_DIR}/lib"
        PDB_OUTPUT_DIRECTORY         "${CMAKE_BINARY_DIR}/lib"
        COMPILE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    set_target_properties(gtest_main
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY     "${CMAKE_BINARY_DIR}/bin"
        LIBRARY_OUTPUT_DIRECTORY     "${CMAKE_BINARY_DIR}/lib"
        ARCHIVE_OUTPUT_DIRECTORY     "${CMAKE_BINARY_DIR}/lib"
        PDB_OUTPUT_DIRECTORY         "${CMAKE_BINARY_DIR}/lib"
        COMPILE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
endif()
