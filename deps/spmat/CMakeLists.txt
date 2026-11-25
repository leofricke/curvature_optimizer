add_library(spmat INTERFACE)
add_library(spmat::spmat ALIAS spmat)

target_include_directories(spmat INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(spmat INTERFACE Eigen3::Eigen)

if(APPLE AND NOT DEFINED OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
    find_program(BREW_EXECUTABLE brew)
    if(BREW_EXECUTABLE)
        execute_process(
            COMMAND ${BREW_EXECUTABLE} --prefix libomp
            OUTPUT_VARIABLE LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE LIBOMP_RESULT
            ERROR_QUIET
        )
        if(LIBOMP_RESULT EQUAL 0 AND EXISTS "${LIBOMP_PREFIX}/include/omp.h")
            set(OpenMP_ROOT "${LIBOMP_PREFIX}")
        endif()
    endif()
endif()

option(SPMAT_REQUIRE_OPENMP "Fail configuration if OpenMP is not found" OFF)
if(SPMAT_REQUIRE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
else()
    find_package(OpenMP COMPONENTS CXX)
endif()
if(OpenMP_CXX_FOUND)
    target_link_libraries(spmat INTERFACE OpenMP::OpenMP_CXX)
endif()
