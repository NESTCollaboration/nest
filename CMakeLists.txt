
# Set minimum CMake version required and prevent building in source
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MACOSX_RPATH 1)
project(NEST)

# Optional GEANT4 integration and ROOT compilation
option(G4 "Build integration with Geant4" OFF)
option(BUILD_ROOT "Build ROOT tools" OFF)

# If the install directory not specified, use the build directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   MESSAGE(STATUS "Using default installation directory in build directory")
 set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif ()

# Set the directory for CMake files, header files, and libraries
set(INSTALL_LIB_DIR "lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR "include" CACHE PATH
  "Installation directory for header files")
set(INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH
  "Installation directory for CMake files")

# Properly set relative paths to the install directory
foreach(p LIB INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Check that the specified compiler has C++11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    add_compile_options(-std=c++11)
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Include all directories which contain needed header files
include_directories("${PROJECT_SOURCE_DIR}")  
include_directories("${PROJECT_SOURCE_DIR}/include")  
include_directories("${PROJECT_SOURCE_DIR}/Detectors")  

# Gather implementation and header files
file(GLOB core ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/Detectors/*.hh)

# Create NEST as a static library with specified headers
add_library(NEST_core STATIC ${core} ${headers} )

# Add a testNEST executable and link to NEST library
add_executable(testNEST testNEST.cpp)
target_link_libraries(testNEST NEST_core)

# Add a bareNEST executable and link to NEST library
add_executable(bareNEST ${PROJECT_SOURCE_DIR}/Tools/bareNEST.cpp)
target_link_libraries(bareNEST NEST_core)

# If doing GEANT4 integration, link GEANT4 package to the NEST library
if(G4)
    find_package(Geant4 REQUIRED)
    include_directories(${Geant4_INCLUDE_DIR})
    file(GLOB g4src ${PROJECT_SOURCE_DIR}/G4integration/*.cc)
    file(GLOB g4headers ${PROJECT_SOURCE_DIR}/G4integration/*.hh)
    add_library(NEST SHARED ${core} ${g4src} ${headers} ${g4headers})
    set_target_properties(NEST PROPERTIES PUBLIC_HEADER "${headers};${g4headers}")
else(G4)
    add_library(NEST SHARED ${core} ${headers})
    set_target_properties(NEST PROPERTIES PUBLIC_HEADER "${headers}")
endif()

# If doing ROOT compilation, add a rootNEST executable (install the various executables either way)
if(BUILD_ROOT)
	list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} $ENV{ROOTSYS}/etc/cmake $ENV{ROOTSYS}/etc/root/cmake)
	list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS} $ENV{ROOTSYS}/etc/cmake $ENV{ROOTSYS}/etc/root/cmake)
	find_package(ROOT REQUIRED)
	include_directories(${ROOT_INCLUDE_DIR})
	include_directories(${ROOT_INCLUDE_DIRS})
	add_executable(rootNEST ${PROJECT_SOURCE_DIR}/Tools/rootNEST.cpp)
	target_link_libraries(rootNEST NEST_core ${ROOT_LIBRARIES})
	install(TARGETS NEST testNEST bareNEST rootNEST
		# IMPORTANT: Add the library to the "export-set"
		EXPORT NESTTargets
					RUNTIME DESTINATION "${INSTALL_CMAKE_DIR}"
		LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
		PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")
else(BUILD_ROOT)
	install(TARGETS NEST testNEST bareNEST
		# IMPORTANT: Add the library to the "export-set"
		EXPORT NESTTargets
					RUNTIME DESTINATION "${INSTALL_CMAKE_DIR}"
		LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
		PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")
endif()



# Define the NEST package and populate CMake files
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
