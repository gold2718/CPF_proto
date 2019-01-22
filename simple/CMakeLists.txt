# Set default project to unknown
if(NOT PROJECT)
  message(STATUS "Setting CCPP project to 'unknown' as none was specified.")
  set(PROJECT "Unknown")
endif (NOT PROJECT)

#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.11)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
    project(simple_cam VERSION 1.0.0)
else(POLICY CMP0048)
    project(simple_cam)
    set(PROJECT_VERSION 1.0.0)
    set(PROJECT_VERSION_MAJOR 1)
    set(PROJECT_VERSION_MINOR 0)
    set(PROJECT_VERSION_PATCH 0)
endif(POLICY CMP0048)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

#------------------------------------------------------------------------------
set(PACKAGE "cam_kessler")
set(AUTHORS  "Cheryl Craig" "Steve Goldhaber" )

#------------------------------------------------------------------------------
# Enable Fortran
enable_language(Fortran)

# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ccpp-framework/cmake")

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
else(OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif()

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# By default we want a static library
option(BUILD_SHARED_LIBS "Build a shared library" OFF)

#------------------------------------------------------------------------------
# Add the CCPP include/module directory
set(CCPP_INCLUDE_DIRS "" CACHE FILEPATH "Path to ccpp includes")
set_property(DIRECTORY PROPERTY INCLUDE_DIRECTORIES ${CCPP_INCLUDE_DIRS})

#------------------------------------------------------------------------------
# Add the CCPP library
set(CCPP_LIB_DIRS "" CACHE FILEPATH "Path to ccpp library")
link_directories(${CCPP_LIB_DIRS})
list(APPEND LIBS "ccpp")

#------------------------------------------------------------------------------
# Set the sources: physics schemes
set(BASE_FILES
  ./src/machine.F
  )
set(SCHEMES
  ./src/kessler.F90
  ./src/converters.F90
#      ./src/kessler_update.F90
  )

# Read the list of files from capfiles.txt and set to CAPS
file(STRINGS ${CCPP_CAPS}/capfiles.txt CAPS)
  message ("CCPP_CAPS=" ${CCPP_CAPS})

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${GFS_PHYS_SRC}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

#apply general fortran tags to all fortran source files
if(${CMAKE_VERSION} LESS 3.3)
  string (REPLACE ";" " " f_flags_str "${f_flags}")
  SET_PROPERTY(SOURCE ${BASE_FILES} ${SCHEMES} ${CAPS} APPEND_STRING PROPERTY COMPILE_FLAGS " ${f_flags_str}")
else(${CMAKE_VERSION} LESS 3.3)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
endif (${CMAKE_VERSION} LESS 3.3)

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/simple)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/library)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------

add_library(base ${BASE_FILES})
add_library(simple_cam ${SCHEMES} ${CAPS})
target_link_libraries(simple_cam LINK_PUBLIC ${LIBS} )
set_target_properties(base PROPERTIES VERSION ${PROJECT_VERSION}
                           SOVERSION ${PROJECT_VERSION_MAJOR}
                           COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                           LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
set_target_properties(simple_cam PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR}
                                     COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                                     LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
# DH* hack for PGI compiler: rename objects in scheme cap object files for ISO_C compliancy
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
set(CAPOBJS)
foreach(cap ${CAPS})
    string(REPLACE "_cap.F90" "_cap.F90.o" capobj "./${CMAKE_FILES_DIRECTORY}/simple_cam.dir/${cap}")
    list(APPEND CAPOBJS ${capobj})
endforeach(cap)

add_custom_command(TARGET simple_cam
                   PRE_LINK
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/pgifix.py --cmake ${CAPOBJS}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Running pgifix_wrapper.py over all scheme caps")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
# *DH end hack for PGI compiler
