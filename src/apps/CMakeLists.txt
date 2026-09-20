function(add_app name sources)
    add_executable(${name} ${sources})

    if(EMSCRIPTEN)
        set_target_properties(${name} PROPERTIES
            LINK_FLAGS "--shell-file ${PROJECT_SOURCE_DIR}/src/viewer/data/shell.html"
            LINK_DEPENDS ${PROJECT_SOURCE_DIR}/src/viewer/data/shell.html
        )
    endif()
endfunction()

add_subdirectory(curvature_optimizer)