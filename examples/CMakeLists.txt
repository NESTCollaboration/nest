
set(MIN_ROOT_VERSION 6.16)


# bareNEST

add_executable(bareNEST bareNEST.cpp)
target_link_libraries(bareNEST PUBLIC NESTCore)

set(EXAMPLES bareNEST)

add_subdirectory("LArNEST/")

# ROOT example
if(BUILD_ROOT)
# https://cliutils.gitlab.io/modern-cmake/examples/root-simple/
find_package(ROOT ${MIN_ROOT_VERSION} REQUIRED CONFIG)

add_executable(rootNEST ${PROJECT_SOURCE_DIR}/examples/rootNEST.cpp)

# List of ROOT targets: https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
# target_link_libraries(rootNEST NEST::Core ROOT::Hist ROOT::Core)
# above fails to find libtbb

include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)
target_link_libraries(rootNEST PUBLIC NESTCore ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})


list(APPEND EXAMPLES rootNEST)
endif()

add_executable(pulseShape pulseShape.cpp)
target_link_libraries(pulseShape PUBLIC NESTCore)
list(APPEND EXAMPLES pulseShape)

add_executable(multipleScatter multipleScatter.cpp)
target_link_libraries(multipleScatter PUBLIC NESTCore)
list(APPEND EXAMPLES multipleScatter)

install(
    TARGETS ${EXAMPLES}        
    RUNTIME DESTINATION bin
)
