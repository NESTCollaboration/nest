
# Set minimum CMake version required and prevent building in source
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MACOSX_RPATH 1)
cmake_policy(SET CMP0048 NEW)
project(NEST2 VERSION 2.0.1)

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
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include all directories which contain needed header files
#target_include_directories("${PROJECT_SOURCE_DIR}/include")
#
#target_include_directories("${PROJECT_SOURCE_DIR}/include/Detectors")

# Gather implementation and header files
file(GLOB core ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/NEST/*.hh ${PROJECT_SOURCE_DIR}/include/Detectors/*.hh)

# Create NEST as a static library with specified headers
add_library(NEST_core STATIC ${core} )
target_include_directories(NEST_core PUBLIC ${PROJECT_SOURCE_DIR}/include/NEST/ ${PROJECT_SOURCE_DIR}/include/Detectors/)

# Add a testNEST executable and link to NEST library
add_executable(testNEST ${PROJECT_SOURCE_DIR}/src/testNEST.cpp)
target_link_libraries(testNEST NEST_core)

# Add a bareNEST executable and link to NEST library
add_executable(bareNEST ${PROJECT_SOURCE_DIR}/examples/bareNEST.cpp)
target_link_libraries(bareNEST NEST_core)

# If doing GEANT4 integration, link GEANT4 package to the NEST library
if(G4)
    find_package(Geant4 REQUIRED)
    include_directories(${Geant4_INCLUDE_DIR})
    file(GLOB g4src ${PROJECT_SOURCE_DIR}/G4integration/*.cpp)
    add_library(NEST SHARED ${core} ${g4src})
    target_link_libraries(NEST ${Geant4_LIBRARIES})
    target_include_directories(NEST PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/G4integration/>
                                           $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>)
    file(GLOB g4headers ${PROJECT_SOURCE_DIR}/G4integration/*.hh)
    list(APPEND headers ${g4headers})
else(G4)
    add_library(NEST SHARED ${core})
endif()
target_include_directories(NEST PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/NEST/;${PROJECT_SOURCE_DIR}/include/Detectors/>"
                                       $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>                 )
set_target_properties(NEST PROPERTIES PUBLIC_HEADER "${headers}")


# If doing ROOT compilation, add a rootNEST executable (install the various executables either way)
if(BUILD_ROOT)
	list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} $ENV{ROOTSYS}/etc/cmake $ENV{ROOTSYS}/etc/root/cmake)
	list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS} $ENV{ROOTSYS}/etc/cmake $ENV{ROOTSYS}/etc/root/cmake)
	find_package(ROOT REQUIRED)
	include_directories(${ROOT_INCLUDE_DIR})
	include_directories(${ROOT_INCLUDE_DIRS})
	add_executable(rootNEST ${PROJECT_SOURCE_DIR}/examples/rootNEST.cpp)
	target_link_libraries(rootNEST NEST_core ${ROOT_LIBRARIES})
	install(TARGETS rootNEST
		RUNTIME DESTINATION "${INSTALL_CMAKE_DIR}"
		LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
		PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")
endif()
install(TARGETS NEST testNEST bareNEST
        # IMPORTANT: Add the library to the "export-set"
        EXPORT NESTConfig
        RUNTIME DESTINATION "${INSTALL_CMAKE_DIR}"
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
        PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")





# Define the NEST package and populate CMake files
export(PACKAGE NEST)
export(TARGETS NEST
  FILE "${PROJECT_BINARY_DIR}/NESTConfig.cmake")


include(CMakePackageConfigHelpers)

write_basic_package_version_file("${PROJECT_BINARY_DIR}/NESTConfigVersion.cmake"
                                 COMPATIBILITY ExactVersion )
install(FILES "${PROJECT_BINARY_DIR}/NESTConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" )

# Install the export set for use with the install-tree
install(EXPORT NESTConfig DESTINATION
  "${INSTALL_CMAKE_DIR}")
