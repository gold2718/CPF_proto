CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(cam_kessler)
ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

# Set the default CAM model to be CAM7
SET(CAM_VERSION "CAM7" CACHE STRING "CAM model version (default: 'CAM7')")
string(TOLOWER "${CAM_VERSION}" CAM_VERSTRING)
# By default, no verbose output
SET(VERBOSITY 0 CACHE STRING "Verbosity level of output (default: 0)")
# By default, generated caps go in ccpp subdir
SET(CCPP_CAP_FILES "${CMAKE_BINARY_DIR}/ccpp" CACHE
  STRING "Location of CCPP-generated cap files")

#MESSAGE("The NETCDF env. variable is $ENV{NETCDF}")
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{NETCDF})

get_filename_component(CAM_DRIVER "${CMAKE_SOURCE_DIR}" DIRECTORY)
get_filename_component(CPF_ROOT ${CAM_DRIVER} DIRECTORY)
SET(CAM_SRC        ${CAM_DRIVER}/src)
SET(SUITES         ${CAM_DRIVER}/suites)
SET(CCPP_FRAMEWORK ${CPF_ROOT}/ccpp-framework)
SET(SIMPLE         ${CPF_ROOT}/simple)
SET(SIMPLE_SRC     ${SIMPLE}/src)
SET(MICM_CHEM      ${CPF_ROOT}/MICM_chemistry)
SET(MICM_CHEM_SRC  ${MICM_CHEM}/src)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    #message(STATUS "Setting build type to 'Debug' as none was specified.")
    #set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                 "MinSizeRel" "RelWithDebInfo")
endif()

ADD_COMPILE_OPTIONS(-O0)

SET(kessler_source_files "${CAM_VERSTRING}_kessler.F90")

SET(library_files
#            shr_infnan_mod.F90
            physconst.F90
            physics_types_${CAM_VERSTRING}.F90
            constituents.F90
            ppgrid.F90
            time_manager.F90
)

#this should call the custome FindNetCDF.cmake file in /src/cmake/modules and populate the NETCDF_INCLUDES and NETCDF_LIBRARIES variables with the appropriate values
SET(NETCDF_F90 "YES")
FIND_PACKAGE(NetCDF)

INCLUDE_DIRECTORIES(${NETCDF_INCLUDES})

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
#  gfortran
#  MESSAGE("gfortran being used.")
  ADD_COMPILE_OPTIONS(-fcheck=all)
  ADD_COMPILE_OPTIONS(-fbacktrace)
  ADD_COMPILE_OPTIONS(-ffpe-trap=zero)
  ADD_COMPILE_OPTIONS(-finit-real=nan)
  ADD_COMPILE_OPTIONS(-ggdb)
  ADD_COMPILE_OPTIONS(-ffree-line-length-none)
  ADD_COMPILE_OPTIONS(-cpp)
  #set(f_flags -fcheck=all -fbacktrace -ffpe-trap=zero -finit-real=nan -ggdb -ffree-line-length-none -g -cpp)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
#  ifort
#  MESSAGE("ifort being used.")
  #ADD_COMPILE_OPTIONS(-check all)
  ADD_COMPILE_OPTIONS(-fpe0)
  ADD_COMPILE_OPTIONS(-warn)
  ADD_COMPILE_OPTIONS(-traceback)
  ADD_COMPILE_OPTIONS(-debug extended)
  ADD_COMPILE_OPTIONS(-fpp)
  #set(f_flags -g -fpe0 -warn -traceback -debug extended -fpp)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
#  pgf90
#  MESSAGE("pgf90 being used.")
  ADD_COMPILE_OPTIONS(-g)
  ADD_COMPILE_OPTIONS(-Mipa=noconst)
  ADD_COMPILE_OPTIONS(-traceback)
  ADD_COMPILE_OPTIONS(-Mfree)
  ADD_COMPILE_OPTIONS(-Mfptrap)
  ADD_COMPILE_OPTIONS(-Mpreprocess)
  #set(f_flags -g -Mipa=noconst -traceback -Mfree -Mfptrap -Mpreprocess)
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#------------------------------------------------------------------------------
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../ccpp-framework/cmake")
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

add_library(library ${library_files})

file(READ "${SIMPLE_SRC}/simple_schemes.txt" FILE_LIST)
file(WRITE "${SIMPLE_SRC}/all_schemes.txt" ${FILE_LIST})
file(READ "${MICM_CHEM_SRC}/MICM_chemistry_schemes.txt" FILE_LIST)
file(APPEND "${SIMPLE_SRC}/all_schemes.txt" ${FILE_LIST})

# Run ccpp_capgen
set(CAPGEN_CMD "${CCPP_FRAMEWORK}/scripts/ccpp_capgen.py")
set(MODEL_FILES "${CAM_VERSTRING}_ccpp_files.txt")
list(APPEND CAPGEN_CMD "${CAM_SRC}/${MODEL_FILES}")
list(APPEND CAPGEN_CMD "${SIMPLE_SRC}/all_schemes.txt")
list(APPEND CAPGEN_CMD "--sdf-pathlist")
#list(APPEND CAPGEN_CMD "${SUITES}/suite_cam_kessler_test_simple1.xml")
list(APPEND CAPGEN_CMD "${SUITES}/suite_kessler_micm.xml")
list(APPEND CAPGEN_CMD "--output-root")
list(APPEND CAPGEN_CMD "${CCPP_CAP_FILES}")
list(APPEND CAPGEN_CMD "--generate-host-cap")
while (VERBOSITY GREATER 0)
  list(APPEND CAPGEN_CMD "--verbose")
  MATH(EXPR VERBOSITY "${VERBOSITY} - 1")
endwhile ()
#MESSAGE(STATUS "capgen cmd: ${CAPGEN_CMD}")
EXECUTE_PROCESS(COMMAND ${CAPGEN_CMD} WORKING_DIRECTORY ${CPF_ROOT}
  OUTPUT_VARIABLE CAPGEN_OUT ERROR_VARIABLE CAPGEN_OUT RESULT_VARIABLE RES)
MESSAGE(STATUS "${CAPGEN_OUT}")
if (RES EQUAL 0)
  MESSAGE(STATUS "CCPP cap generation completed")
else(RES EQUAL 0)
  MESSAGE(FATAL_ERROR "CCPP cap generation FAILED: result = ${RES}")
endif(RES EQUAL 0)

# Read the list of files from capfiles.txt and set to CCPP_CAPS
file(STRINGS ${CCPP_CAP_FILES}/capfiles.txt CCPP_CAPS)
add_library(caps ${CCPP_CAPS})
INCLUDE_DIRECTORIES(${CCPP_CAP_FILES})

ADD_SUBDIRECTORY(${SIMPLE} ${CMAKE_BINARY_DIR}/simple)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/simple)

ADD_SUBDIRECTORY(${MICM_CHEM} ${CMAKE_BINARY_DIR}/MICM_chemistry)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/MICM_chemistry)

ADD_DEPENDENCIES(library base)
ADD_DEPENDENCIES(simple_cam library)
ADD_DEPENDENCIES(micmchem library)
ADD_DEPENDENCIES(caps simple_cam)
ADD_DEPENDENCIES(caps micmchem)

ADD_EXECUTABLE(cam_kessler ${kessler_source_files})
TARGET_LINK_LIBRARIES(cam_kessler caps simple_cam micmchem library base ${NETCDF_LIBRARIES})
set_target_properties(cam_kessler PROPERTIES
                               COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                               LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
