if(NOT TARGET gcem)
    include(FetchContent)

    FetchContent_Declare(
            gcem
            GIT_REPOSITORY https://github.com/kthohr/gcem
            GIT_TAG        v1.13.1
    )


    FetchContent_GetProperties(gcem)
    if(NOT gcem_POPULATED)
        FetchContent_Populate(gcem)
        add_subdirectory(${gcem_SOURCE_DIR} ${gcem_BINARY_DIR})
    endif()
endif()
