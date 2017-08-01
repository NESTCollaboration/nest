
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
project(NEST)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   MESSAGE(STATUS "Using default installation directory in build directory")
 set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif ()

set(INSTALL_LIB_DIR "lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR "include" CACHE PATH
  "Installation directory for header files")
set(INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH
  "Installation directory for CMake files")

foreach(p LIB INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    add_compile_options(-std=c++11)
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

include_directories("${PROJECT_SOURCE_DIR}")  


add_library(NEST SHARED NEST.cpp NEST.hh)
set_target_properties(NEST PROPERTIES PUBLIC_HEADER "NEST.hh")
install(TARGETS NEST
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT NESTTargets
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")

add_executable(testNEST testNEST.cpp)
target_link_libraries(testNEST NEST)

export(PACKAGE NEST)
export(TARGETS NEST
  FILE "${PROJECT_BINARY_DIR}/NESTTargets.cmake")


file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(NESTConfig.cmake.in
  "${PROJECT_BINARY_DIR}/NESTConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${NEST_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(NESTConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/NESTConfig.cmake" @ONLY)
# ... for both
configure_file(NESTConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/NESTConfigVersion.cmake" @ONLY)
 
# Install the NESTConfig.cmake and NESTConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/NESTConfig.cmake"
  "${PROJECT_BINARY_DIR}/NESTConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT NESTTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)