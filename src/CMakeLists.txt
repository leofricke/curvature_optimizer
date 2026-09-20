if (WIN32)
    add_compile_options(/W4)
    add_compile_options(/wd4267) # conversion from size_t
    add_compile_options(/wd4244) # conversion from double to float
    add_compile_options(/wd4305) # truncation from double to float

    if(MSVC)
        add_compile_options(-D_USE_MATH_DEFINES -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS)
    else()
        add_compile_options(-D_USE_MATH_DEFINES -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS -fopenmp)
    endif()
else ()
    add_compile_options(-Wall)
    add_compile_options(-Wextra)
    add_compile_options(-Wpedantic)
    add_compile_options(-Wshadow)
endif()

option(ENABLE_ASAN "Build with AddressSanitizer" OFF)
if(ENABLE_ASAN AND UNIX)
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
    message(STATUS "AddressSanitizer is enabled")
endif()

if(EMSCRIPTEN)
    add_compile_options(-fexceptions)
    add_compile_options(-pthread)
    add_link_options(
        "SHELL:-sWASM=1 -pthread -sUSE_WEBGL2=1 -sUSE_GLFW=3 -sALLOW_MEMORY_GROWTH=1 -sALLOW_TABLE_GROWTH=1 -sSTACK_SIZE=5MB -sEXPORTED_RUNTIME_METHODS=[ccall] -sNO_DISABLE_EXCEPTION_CATCHING"
    )
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
endif()

add_subdirectory(utils)
add_subdirectory(io)
add_subdirectory(opt)
add_subdirectory(curvature_optimizer)

if(NOT EMSCRIPTEN AND BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(NOT BUILD_HEADLESS)
    add_subdirectory(viewer)
endif()

add_subdirectory(apps)