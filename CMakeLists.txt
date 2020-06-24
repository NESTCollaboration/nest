
# Set minimum CMake version required and prevent building in source
cmake_minimum_required(VERSION 3.8...3.16 FATAL_ERROR)

project(NEST VERSION 2.1.0 LANGUAGES CXX)

option(G4 "Build integration with Geant4" OFF)
option(BUILD_ROOT "Build ROOT tools" OFF)
option(BUILD_EXAMPLES "Build NEST examples" ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)


set(NEST_CORE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NEST.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RandomGen.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestSpectra.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/VDetector.cpp
)
set(NEST_CORE_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/NEST.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/RandomGen.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/TestSpectra.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/Detectors/VDetector.hh
)

add_library(Core ${NEST_CORE_SOURCES} ${NEST_CORE_HEADERS})
add_library(NEST::Core ALIAS Core)
list(APPEND NEST_TARGETS "Core")
list(APPEND NEST_LIBRARY_TARGETS "NEST::Core")

target_include_directories(
  Core PUBLIC
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/NEST>
  $<INSTALL_INTERFACE:include/Detectors>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/NEST>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Detectors>
)

target_compile_features(Core PUBLIC cxx_std_11)

add_executable(testNEST
  ${CMAKE_CURRENT_SOURCE_DIR}/src/testNEST.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/testNEST.hh
)
target_link_libraries(testNEST PUBLIC NEST::Core)
list(APPEND NEST_TARGETS "testNEST")

## Optional targets

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(G4)
    add_subdirectory(G4integration)
    list(APPEND NEST_TARGETS "NESTG4")
    set(NEST_LIBRARY_TARGETS "${NEST_LIBRARY_TARGETS} NEST::NESTG4")
endif()

# installation, Targets.cmake and Config.cmake
message(STATUS "NEST_LIBRARY_TARGETS: ${NEST_LIBRARY_TARGETS}")
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/NEST)

write_basic_package_version_file(
    NESTConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY ExactVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/NESTConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/NESTConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/NESTConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/NESTConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_BINDIR}")
message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
install(
  TARGETS ${NEST_TARGETS}
  EXPORT NEST-Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT NEST-Targets
  FILE NESTTargets.cmake
  NAMESPACE NEST::
  DESTINATION ${INSTALL_CONFIGDIR}
)

##############################################

export(
  EXPORT NEST-Targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/NESTTargets.cmake
  NAMESPACE NEST::
)

# Register package in user's package registry
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE NEST)
