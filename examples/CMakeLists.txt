
set(MIN_ROOT_VERSION 6.16)


# bareNEST

add_executable(bareNEST bareNEST.cpp)
target_link_libraries(bareNEST PUBLIC NEST::Core)

set(EXAMPLES bareNEST)

# ROOT example
if(BUILD_ROOT)
# https://cliutils.gitlab.io/modern-cmake/examples/root-simple/
find_package(ROOT ${MIN_BOOST_VERSION} REQUIRED CONFIG)

add_executable(rootNEST ${PROJECT_SOURCE_DIR}/examples/rootNEST.cpp)

# List of ROOT targets: https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
# target_link_libraries(rootNEST NEST::Core ROOT::Hist ROOT::Core)
# above fails to find libtbb

include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)
target_link_libraries(rootNEST PUBLIC NEST::Core ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})


list(APPEND EXAMPLES rootNEST)
endif()


install(
    TARGETS ${EXAMPLES}        
    RUNTIME DESTINATION bin
)
