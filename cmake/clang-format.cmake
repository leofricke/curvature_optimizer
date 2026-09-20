# find clang-format executable
find_program(
    CLANG_FORMAT_EXE
    NAMES "clang-format"
    DOC "Path to clang-format executable")

# add custom target if found
if(CLANG_FORMAT_EXE)
    message(STATUS "clang-format found: ${CLANG_FORMAT_EXE}")

    file(
        GLOB_RECURSE
        ALL_SOURCES
        "src/**/*.cpp"
        "src/**/*.h"
    )

    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXE}
        -style=file
        -i ${ALL_SOURCES}
        COMMENT "Running clang-format"
    )
else()
    message(STATUS "clang-format not found.")
endif()