
# Set minimum CMake version required and prevent building in source
cmake_minimum_required(VERSION 3.8...3.16 FATAL_ERROR)

###############        Get version from git        #####################
execute_process(COMMAND git describe --tag
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE PROJECT_GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)

function(get_versions versionString version)
    if ("${versionString}" STREQUAL "")
        set(version "0.0.0" PARENT_SCOPE)
        return()
    endif ()

    string(REGEX REPLACE "v([0-9]*)([.][0-9]*[.][0-9]*-?.*)$" "\\1" numbers ${versionString})
    set(major ${numbers})
    string(REGEX REPLACE "v([0-9]*[.])([0-9]*)([.][0-9]*-?.*)$" "\\2" numbers ${versionString})
    set(minor ${numbers})
    string(REGEX REPLACE "v([0-9]*[.][0-9]*[.])([0-9]*)(-?.*)$" "\\2" numbers ${versionString})
    set(patch ${numbers})
    set(version "${major}.${minor}.${patch}" PARENT_SCOPE)
endfunction()

get_versions("${PROJECT_GIT_VERSION}" version)
set(PROJECT_GIT_VERSION ${version})

###############    project name and version        #####################
project(NEST VERSION ${PROJECT_GIT_VERSION} LANGUAGES CXX)

option(G4 "Build integration with Geant4" OFF)
option(BUILD_ROOT "Build ROOT tools" OFF)
option(BUILD_EXAMPLES "Build NEST examples" ON)
option(BUILD_GARFIELD "Build integration with Garfield++" ON)

# enable position independent code so shared libraries can link against 
# the static version of NEST::Core 
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(cmake/gcem.cmake)



set(NEST_CORE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NEST.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RandomGen.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestSpectra.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/VDetector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/GammaHandler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidityTests.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LArNEST.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LArDetector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestLArSpectra.cpp
)
set(NEST_CORE_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/NEST.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/RandomGen.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/TestSpectra.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/Detectors/VDetector.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/Detectors/LArDetector.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/GammaHandler.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/ValidityTests.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/LArNEST.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/TestLArSpectra.hh
)

add_library(NESTCore ${NEST_CORE_SOURCES} ${NEST_CORE_HEADERS})

target_link_libraries(NESTCore PUBLIC gcem)

list(APPEND NEST_TARGETS "NESTCore")

target_include_directories(
  NESTCore PUBLIC
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/NEST>
  $<INSTALL_INTERFACE:include/Detectors>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/NEST>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Detectors>
)

target_compile_features(NESTCore PUBLIC cxx_std_17)

add_executable(execNEST
  ${CMAKE_CURRENT_SOURCE_DIR}/src/execNEST.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/include/NEST/execNEST.hh
)
target_link_libraries(execNEST PUBLIC NESTCore)
list(APPEND NEST_TARGETS "execNEST")

## scripts
configure_file(examples/loopNEST.in loopNEST)
configure_file(examples/220RnCalib.in 220RnCalib)

## Optional targets

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(BUILD_GARFIELD)
  add_subdirectory(GarfieldppIntegration)
endif()

if(G4)
    add_subdirectory(G4integration)
    list(APPEND NEST_TARGETS "NESTG4")
endif()

# installation, Targets.cmake and Config.cmake
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

## install scripts

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/loopNEST
  ${CMAKE_CURRENT_BINARY_DIR}/220RnCalib
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  DESTINATION ${CMAKE_INSTALL_BINDIR}
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
