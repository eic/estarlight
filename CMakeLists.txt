###########################################################################
#
#    Copyright 2010
#
#    This file is part of Starlight.
#
#    Starlight is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#	  
#    Starlight is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#	  
#    You should have received a copy of the GNU General Public License
#    along with Starlight. If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################
#
# File and Version Information:
# $Rev:: 247                         $: revision of last commit
# $Author:: butter                   $: author of last commit
# $Date:: 2016-03-05 00:59:24 +0000 #$: date of last commit
#
# Description:
#      Starlight build file
#
#
###########################################################################


# check if cmake has the required version
cmake_minimum_required(VERSION 3.6.0 FATAL_ERROR)

# set verbosity
set(CMAKE_VERBOSE_MAKEFILE 0)  # if set to 1 compile and link commands are displayed during build
# the same effect can be achieved by calling 'make VERBOSE=1'


# The version number. 9999 indicates trunk
set (Starlight_VERSION_MAJOR 9999) 
set (Starlight_VERSION_MINOR 1)
set (Starlight_VERSION_MINOR_MINOR 0)

# define project
project(starlight)
#find_package (Threads)

# load some common cmake macros
# set path, where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'")
include(CommonMacros)


# force out-of-source builds.
enforce_out_of_source_build()


# warn user if system is not UNIX
if(NOT UNIX)
  message(FATAL_ERROR "This is an unsupported system.")
endif()
message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version '${CMAKE_HOST_SYSTEM_VERSION}' architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}'")
message(STATUS "Compiling for system '${CMAKE_SYSTEM_NAME}' version '${CMAKE_SYSTEM_VERSION}' architecture '${CMAKE_SYSTEM_PROCESSOR}'")

option (CPP11 "Enable compilation with C++11 features" ON) 

# define build types
# set a default build type for single-configuration CMake generators, if no build type is set.
set(CMAKE_BUILD_TYPE Debug)
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type was specified. Setting build type to 'Release'.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# common compiler flags
if (CMAKE_COMPILER_IS_GNUCC)
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
  message(STATUS "GCC_VERSTION")
  message(STATUS  ${GCC_VERSION})
  if (GCC_VERSION VERSION_GREATER 4.6 OR GCC_VERSION VERSION_EQUAL 4.6)
    message(STATUS "GCC_VERSION>=4.6")
    if(CPP11)
      set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wno-error=unused-but-set-variable -Wno-error=unused-but-set-parameter -std=c++11")
      message(STATUS "Enabling usage of C++11 features")
    else()
      set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wno-error=unused-but-set-variable -Wno-error=unused-but-set-parameter")
    endif()
  else()
    message(STATUS "GCC_VERSION<4.6")
      set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wno-error=unused-but-set-variable -Wno-error=unused-but-set-parameter")
    if(CPP11)
      message(WARNING "C++11 features not supported for your compiler")
    endif()
  endif()
else()
  message(STATUS "Not GCC")
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")
  if(CPP11)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -std=c++11")
    # message(WARNING "C++11 features not supported for your compiler")
    message(WARNING "Trying to enable C++11 features with '-std=c++11'")
  endif()
endif()
# flags for specific build types
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_LDFLAGS_DEBUG "-g")
# report global build settings
message(STATUS "Using CXX compiler '${CMAKE_CXX_COMPILER}'")
message(STATUS "Using CXX general compiler flags '${CMAKE_CXX_FLAGS}'")
foreach(_BUILD_TYPE "DEBUG" "MINSIZEREL" "RELEASE" "RELWITHDEBINFO")
  message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS_${_BUILD_TYPE}}' for build type ${_BUILD_TYPE}")
endforeach()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'")


# redirect output files
#set(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib")
message(STATUS "Using library output path '${LIBRARY_OUTPUT_PATH}'")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")
message(STATUS "Using executable output path '${EXECUTABLE_OUTPUT_PATH}'")


# make CMAKE_SOURCE_DIR accessible in source code via predefined macro CMAKE_SOURCE_DIR
if(CMAKE_SOURCE_DIR)
  add_definitions(-D'CMAKE_SOURCE_DIR=\"${CMAKE_SOURCE_DIR}\"')
else()
  add_definitions(-D'CMAKE_SOURCE_DIR=\"\"')
endif()


# make SVN version string accessible in source code via predefined macro SVN_VERSION
find_package(Subversion)
if(Subversion_FOUND)
  # unfortunately CMAKE only parses 'svn info'
  find_program(SVNVERSION_EXECUTABLE
    svnversion
		)
  if(NOT SVNVERSION_EXECUTABLE)
    message(STATUS "Could not find subversion command 'svnversion'. Repository version unknown.")
  else()
    execute_process(
      COMMAND ${SVNVERSION_EXECUTABLE} "${CMAKE_SOURCE_DIR}"
      OUTPUT_VARIABLE SVN_VERSION
      RESULT_VARIABLE _SVNVERSION_RETURN
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT ${_SVNVERSION_RETURN})
      message(STATUS "Subversion repository revision is '${SVN_VERSION}'")
    else()
      message(STATUS "Error running 'svnversion'. Repository version unknown.")
      set(SVN_VERSION "")
    endif()
  endif()
else()
  message(STATUS "Could not find subversion installation. Repository version unknown.")
endif()
if(SVN_VERSION)
  add_definitions(-D'SVN_VERSION=\"${SVN_VERSION}\"')
else()
  add_definitions(-D'SVN_VERSION=\"\"')
endif()


# setup doxygen
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  message(WARNING "Cannot find Doxygen. No HTML documentation will be generated.")
else()
  set(DOXYGEN_TARGET  "doxygen")
  set(DOXYGEN_DOC_DIR "${CMAKE_SOURCE_DIR}/doxygen")
  set(DOXYGEN_CONF    "${CMAKE_SOURCE_DIR}/starlightDoxyfile.conf")
  message(STATUS "Run 'make ${DOXYGEN_TARGET}' to create Doxygen documentation files in '${DOXYGEN_DOC_DIR}'")
  add_custom_target(${DOXYGEN_TARGET}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF}
    DEPENDS ${DOXYGEN_CONF}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		)
endif()

option (ENABLE_HEPMC3  "Enable compilation against hepmc3 (necessary for hepmc3 output)" OFF) 
if(ENABLE_HEPMC3)
  find_package(HepMC3 HINTS ${HepMC3_DIR} ${HepMC3} )
  if(HepMC3_FOUND)
    message ( "-- HepMC3 library found at ${HEPMC3_LIB} --> Adding HepMC sources and headers" )
    set(optionalLibs ${optionalLibs} ${HEPMC3_LIB})
    option(ENABLE_HEPMC3  "Should we use the HepMC3 library" ON) 
  else()
    message ( WARNING "Hepmc3 not found.")
    option(ENABLE_HEPMC3  "Should we use the HepMC3 library" OFF) 
  endif()
endif()

# setup Pythia 8
option (ENABLE_PYTHIA  "Enable compilation against pythia (necessary for certain processes)" OFF) 
if(ENABLE_PYTHIA)
  find_package(Pythia8)
  if(PYTHIA8_FOUND)
    set(optionalLibs ${optionalLibs} ${PYTHIA8_LIBRARY})
#    find_package(LHAPDF REQUIRED) # implemented for dummy version in Pythia8
#    set(optionalLibs ${optionalLibs} ${LHAPDF_LIBRARIES})#liblhapdfdummy
#    removed from v8.2, if you want to add your own lhapdf set, uncomment
#    find_package + set  - jb 05192015
    option(ENABLE_PYTHIA  "Should we use the Pythia8 library" ON) 
  else() 
    option(ENABLE_PYTHIA  "Should we use the Pythia8 library" OFF) 
  endif()
endif()

# setup Pythia 6
option (ENABLE_PYTHIA6  "Enable compilation against pythia 6 (necessary for certain processes)" OFF) 
if(ENABLE_PYTHIA6)
   find_package(Pythia6 REQUIRED)
   if(PYTHIA6_FOUND)
      set(optionalLibs ${optionalLibs} ${PYTHIA6_LIBRARY})
      option (ENABLE_PYTHIA6  "Enable compilation against pythia 6 (necessary for certain processes)" ON) 
      include_directories(pythia6)
   else(PYTHIA6_FOUND) 
      option (ENABLE_PYTHIA6  "Enable compilation against pythia 6 (necessary for certain processes)" OFF) 
   endif(PYTHIA6_FOUND)
endif()

# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_SOURCE_DIR}/include
	${PROJECT_BINARY_DIR}
	${PYTHIA8_INCLUDE_DIR}#uncommented 05192015
	)
include_directories(${INCLUDE_DIRECTORIES})

# Set our source files, include the generated dictionary
set(SOURCES
  src/bessel.cpp
  src/beam.cpp
  src/inputParameters.cpp
  src/beambeamsystem.cpp
  src/starlightparticle.cpp
  src/gammaaluminosity.cpp
  src/randomgenerator.cpp
  src/nucleus.cpp
  src/eventchannel.cpp
  src/gammaavm.cpp
  src/gammagammasingle.cpp
  src/photonNucleusCrossSection.cpp
  src/wideResonanceCrossSection.cpp
  src/narrowResonanceCrossSection.cpp
  src/readinluminosity.cpp
  src/twophotonluminosity.cpp
  src/vector3.cpp
  src/lorentzvector.cpp
  src/filewriter.cpp
  src/eventfilewriter.cpp
  src/starlightparticlecodes.cpp
  src/nBodyPhaseSpaceGen.cpp
  src/inputParser.cpp
  src/incoherentPhotonNucleusLuminosity.cpp
  src/incoherentVMCrossSection.cpp
# eSTARlight
  src/eXevent.cpp
  src/gammaeluminosity.cpp	
  src/e_wideResonanceCrossSection.cpp
  src/e_narrowResonanceCrossSection.cpp
  src/e_starlight.cpp
  src/e_starlightStandalone.cpp
  )

if (ENABLE_HEPMC3)
  set (SOURCES
		${SOURCES}
                src/hepmc3writer.cpp
		)
  include_directories(${HEPMC3_INCLUDE_DIR})
  add_definitions( -DHEPMC3_ON )
endif()
if(ENABLE_PYTHIA) 
  set (SOURCES
		${SOURCES}
		#src/PythiaStarlight.cpp
                src/pythiadecayer.cpp
		)
  include_directories(${PYTHIA8_INCLUDE_DIR})
endif()
if(ENABLE_PYTHIA6)
  set (SOURCES 
		${SOURCES}
		src/starlightpythia.cpp
		src/spectrum.cpp
		src/spectrumprotonnucleus.cpp
		)
endif()

# add Starlight library to the build system
set(THIS_LIB "Starlib")
add_library(${THIS_LIB} STATIC ${SOURCES})
#make_shared_library("${THIS_LIB}" "${SOURCES}"
#	"${PYTHIA8_LIBRARY}"
#	"${LHAPDF_LIBRARIES}"
#)

if(ENABLE_PYTHIA6)
  enable_language(Fortran)
endif()

# add starlight executable to the build system
add_executable(e_starlight src/e_main.cpp)
target_link_libraries(e_starlight Starlib ${optionalLibs})# ${CMAKE_THREAD_LIBS_INIT}) 

configure_file (
  "${PROJECT_SOURCE_DIR}/starlightconfig.h.in"
  "${PROJECT_BINARY_DIR}/starlightconfig.h"
  )
# Erase xsec values in case changes in code affects the xsec, executed during make process
add_custom_command (TARGET Starlib POST_BUILD COMMAND touch ARGS slight.txt)
add_custom_command (TARGET Starlib POST_BUILD COMMAND cp ARGS slight.txt slight.txt.bak)
add_custom_command (TARGET Starlib POST_BUILD COMMAND echo ARGS '' > slight.txt )

# Installation of targets
install(TARGETS e_starlight ${THIS_LIB}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

message(STATUS "Cmake did not find any errors. run 'make' to build the project.")
message(STATUS "On multi-core machines 'make -j#', where # is the number of parallel jobs, can speedup compilation considerably.")
