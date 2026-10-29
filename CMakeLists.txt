# Overall structure
#
# POLICIES
# PROJECT DEFINITIONS
#    Version
#    Project declaration
#    Minimum dependency versions
#    Directories and search paths
#    Packaging
# OPTIONS
#    Declarations with preprocessing/printing
# BASIC SETTINGS
#    Shared/static libraries
#    Build type
#    [LTO]
#    Githash
# DEPENDENCIES
#    Search and process
#    Check for position-independent code
# TOOLS
#    [Code coverage]
#    Global compilation settings
# SUBDIRECTORIES
# INSTALLATION
#    Licenses
#    Config(-version) files
#    Build-tree targets
# [LOGGING]
#    [Non-default parameters]

cmake_minimum_required (VERSION 3.11.0)

##### POLICIES #####

# FindBoost is removed in version cmake 3.30.
if(POLICY CMP0167)
   cmake_policy(SET CMP0167 NEW)
endif()

##### PROJECT DEFINITIONS #####

# PaPILO version
set(PAPILO_VERSION_MAJOR 4)
set(PAPILO_VERSION_MINOR 0)
set(PAPILO_VERSION_PATCH 0)
set(PAPILO_API_VERSION 13)

# Project
project(papilo
   VERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}.${PAPILO_VERSION_PATCH}.0)

set(BOOST_MIN_VERSION 1.67)
set(TBB_MIN_VERSION 2018)

# Set installation directories according to GNU coding standards.
include(GNUInstallDirs)

# Path to CMake modules such as FindGMP.cmake.
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

##### OPTIONS #####

# Options
option(BUILD_PACKAGE "are we building a package" OFF)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

option(GITHASH_BRANCH "append the current branch to the githash" OFF)

set(GMP AUTO CACHE STRING "should PaPILO link to GMP")
set_property(CACHE GMP PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${GMP}" GMP)

set(QUADMATH AUTO CACHE STRING "should quadmath library be linked")
set_property(CACHE QUADMATH PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${QUADMATH}" QUADMATH)

set(LUSOL AUTO CACHE STRING "should LUSOL package be enabled")
set_property(CACHE LUSOL PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${LUSOL}" LUSOL)


set(TBB AUTO CACHE STRING "should TBB be linked")
set_property(CACHE TBB PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${TBB}" TBB)

set(SOPLEX AUTO CACHE STRING "should executable link to the SoPlex LP solver")
set_property(CACHE SOPLEX PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${SOPLEX}" SOPLEX)

set(SYSTEM_FMT OFF CACHE STRING "use system FMT library instead of bundled headers")
set_property(CACHE SYSTEM_FMT PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${SYSTEM_FMT}" SYSTEM_FMT)

set(SCIP AUTO CACHE STRING "should the executable link to the SCIP solver")
set_property(CACHE SCIP PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${SCIP}" SCIP)

set(HIGHS AUTO CACHE STRING "should the executable link to the HiGHS solver")
set_property(CACHE HIGHS PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${HIGHS}" HIGHS)

set(GLOP AUTO CACHE STRING "should the executable link to the GLOP solver from OR-Tools")
set_property(CACHE GLOP PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${GLOP}" GLOP)

set(GUROBI AUTO CACHE STRING "should executable link to the Gurobi solver")
set_property(CACHE GUROBI PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${GUROBI}" GUROBI)

option(BUILD_TESTING "should tests be enabled and built" ON)

option(HEADER_ONLY "turn off compiling specialized templates into library" OFF)

##### BASIC SETTINGS #####

if(BUILD_SHARED_LIBS)
   message(STATUS "Building shared libraries.")
   set(PAPILO_HAS_PIC 1)
else()
   message(STATUS "Building static libraries.")
   set(PAPILO_HAS_PIC 0)
endif()

# Make 'Release' the default build type.
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
elseif(NOT CMAKE_BUILD_TYPE MATCHES "Release|Debug|RelWithDebInfo|MinSizeRel")
   message(FATAL_ERROR "Option CMAKE_BUILD_TYPE has invalid value \"${CMAKE_BUILD_TYPE}\"! Possible are `Release`, `Debug`, `RelWithDebInfo` and `MinSizeRel`.")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Query git hash.
include(GitHash)
githash_query(PAPILO_GITHASH "${CMAKE_CURRENT_SOURCE_DIR}/src/papilo/Config.hpp" "PAPILO_GITHASH" ${GITHASH_BRANCH})
if(PAPILO_GITHASH)
   message(STATUS "Git hash: ${PAPILO_GITHASH}")
   set(PAPILO_HAS_GITHASH 1)
else()
   set(PAPILO_HAS_GITHASH 0)
endif()

# Definition for header-only library.
set(PAPILO_HEADER_ONLY ${HEADER_ONLY})
if(PAPILO_HEADER_ONLY)
   message(STATUS "Creating header-only PaPILO library.")
endif()

##### DEPENDENCIES #####

# Search for system fmt.
set(PAPILO_HAS_SYSTEM_FMT 0)
if(SYSTEM_FMT)
   if(SYSTEM_FMT STREQUAL "AUTO")
      find_package(fmt CONFIG)
   else()
      find_package(fmt CONFIG REQUIRED)
   endif()
   if(fmt_FOUND)
      set(PAPILO_HAS_SYSTEM_FMT 1)
   endif()
endif()
if(PAPILO_HAS_SYSTEM_FMT)
   message(STATUS "Using system FMT.")
else()
   message(STATUS "Using bundled FMT.")
endif()

# Search for Boost: headers are required and the other components are needed for the binaries.
set(PAPILO_HAS_BOOST_IOSTREAMS 0)
set(PAPILO_HAS_BOOST_IOSTREAMS_ZLIB 0)
set(PAPILO_HAS_BOOST_IOSTREAMS_BZIP2 0)
set(PAPILO_HAS_BOOST_PROGRAM_OPTIONS 0)
set(PAPILO_HAS_BOOST_SERIALIZATION 0)
message(STATUS "Searching for Boost.")
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED OPTIONAL_COMPONENTS iostreams program_options serialization)
if(NOT DEFINED Boost_VERSION_STRING)
   math(EXPR _boost_major "${Boost_VERSION} / 100000")
   math(EXPR _boost_minor "${Boost_VERSION} / 100 % 1000")
   math(EXPR _boost_patch "${Boost_VERSION} % 100")
   set(Boost_VERSION_STRING "${_boost_major}.${_boost_minor}.${_boost_patch}")
endif()
message(STATUS "Found Boost version ${Boost_VERSION_STRING}.")
set(PAPILO_HAS_BOOST 1)
set(PAPILO_BOOST_COMPONENTS)

include(CheckCXXSourceCompiles)
if(Boost_iostreams_FOUND)
   set(PAPILO_HAS_BOOST_IOSTREAMS 1)
   set(PAPILO_BOOST_COMPONENTS ${PAPILO_BOOST_COMPONENTS} iostreams)

   set(CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIRS})
   set(CMAKE_REQUIRED_LIBRARIES ${Boost_LIBRARIES})
   set(CMAKE_REQUIRED_QUIET TRUE)
   unset(PAPILO_HAS_BOOST_IOSTREAMS_ZLIB)
   check_cxx_source_compiles(
      "#include <boost/iostreams/filter/gzip.hpp>
       int main() { auto decomp = boost::iostreams::gzip_decompressor(); (void)decomp; return 0; }"
      PAPILO_HAS_BOOST_IOSTREAMS_ZLIB )
   unset(PAPILO_HAS_BOOST_IOSTREAMS_BZIP2)
   check_cxx_source_compiles(
      "#include <boost/iostreams/filter/bzip2.hpp>
       int main() { auto decomp = boost::iostreams::bzip2_decompressor(); (void)decomp; return 0; }"
      PAPILO_HAS_BOOST_IOSTREAMS_BZIP2 )
   unset(CMAKE_REQUIRED_INCLUDES)
   unset(CMAKE_REQUIRED_LIBRARIES)
   unset(CMAKE_REQUIRED_QUIET)
   if(PAPILO_HAS_BOOST_IOSTREAMS_ZLIB AND PAPILO_HAS_BOOST_IOSTREAMS_BZIP2)
      message(STATUS "Tested Boost: support for zlib and bzip2 (de)compression.")
   elseif(PAPILO_HAS_BOOST_IOSTREAMS_ZLIB)
      message(STATUS "Tested Boost: support for zlib (de)compression.")
   elseif(PAPILO_HAS_BOOST_IOSTREAMS_BZIP2)
      message(STATUS "Tested Boost: support for bzip2 (de)compression.")
   else()
      message(STATUS "Tested Boost: no support for zlib or bzip2 (de)compression.")
   endif()
endif()
if(Boost_program_options_FOUND)
   set(PAPILO_HAS_BOOST_PROGRAM_OPTIONS 1)
   set(PAPILO_BOOST_COMPONENTS ${PAPILO_BOOST_COMPONENTS} program_options)
endif()
if(Boost_serialization_FOUND)
   set(PAPILO_HAS_BOOST_SERIALIZATION 1)
   set(PAPILO_BOOST_COMPONENTS ${PAPILO_BOOST_COMPONENTS} serialization)
endif()
if(NOT PAPILO_HAS_BOOST_IOSTREAMS OR NOT PAPILO_HAS_BOOST_PROGRAM_OPTIONS OR NOT PAPILO_HAS_BOOST_SERIALIZATION)
   message(WARNING "Boost iostreams, program_options or serialization not found. PaPILO can only be used as a library!")
endif()

# Search for GMP.
if(GMP)
   if(GMP STREQUAL "AUTO")
      find_package(GMP)
   else()
      find_package(GMP REQUIRED)
   endif()
   if(GMP_FOUND)
      set(PAPILO_HAS_GMP 1)
   else()
      set(PAPILO_HAS_GMP 0)
   endif()
else()
   message(STATUS "GMP is disabled.")
   set(PAPILO_HAS_GMP 0)
endif()

# Search for Quadmath.
if(QUADMATH)
   if(QUADMATH STREQUAL "AUTO")
      find_package(Quadmath)
   else()
      find_package(Quadmath REQUIRED)
   endif()
   if(Quadmath_FOUND)
      set(PAPILO_HAS_QUADMATH 1)
   else()
      set(PAPILO_HAS_QUADMATH 0)
   endif()
else()
   message(STATUS "Quadmath is disabled.")
   set(PAPILO_HAS_QUADMATH 0)
endif()

# Search for TBB.
if(TBB)
   if(TBB STREQUAL "AUTO")
      find_package(TBB ${TBB_MIN_VERSION} HINTS "${TBB_DIR}" COMPONENTS tbb tbbmalloc)
      find_package(Threads)
   else()
      find_package(TBB ${TBB_MIN_VERSION} HINTS "${TBB_DIR}" REQUIRED COMPONENTS tbb tbbmalloc)
      find_package(Threads REQUIRED)
   endif()
   if(TBB_FOUND AND Threads_FOUND)
      message(STATUS "Found TBB version ${TBB_VERSION}: ${TBB_CONFIG}")
      set(PAPILO_HAS_TBB 1)
   elseif(NOT TBB_FOUND)
      message(WARNING "TBB not found. Disabling parallel presolving capabilities!")
      set(PAPILO_HAS_TBB 0)
   else()
      message(WARNING "Threads library not found, but required for TBB. Disabling parallel presolving capabilities!")
      set(PAPILO_HAS_TBB 0)
   endif()
else()
   set(PAPILO_HAS_TBB 0)
   message(WARNING "Parallel presolving capabilities disabled due to TBB being disabled!")
endif()

# LUSOL
if(LUSOL)
   include(CheckLanguage)
   check_language(Fortran)
   if(LUSOL STREQUAL "AUTO")
      if(CMAKE_Fortran_COMPILER)
         find_package(BLAS)
      else()
         message(STATUS "No Fortran support.")
      endif()
   elseif(NOT CMAKE_Fortran_COMPILER)
      message(FATAL_ERROR "Fortran support is required due to LUSOL=${LUSOL}.")
   else()
      find_package(BLAS REQUIRED)
   endif()
   if(CMAKE_Fortran_COMPILER AND BLAS_FOUND)
      enable_language(Fortran)
      set(PAPILO_HAS_LUSOL 1)
      if(CMAKE_VERSION VERSION_LESS 3.18)
         # TODO: Remove this once we require CMake >= 3.18.
         message(WARNING "Defining target BLAS::BLAS with linker flags \"${BLAS_LINKER_FLAGS}\" and libraries \"${BLAS_LIBRARIES}\". Consider upgrading to CMake version 3.18 or higher.")
         add_library(BLAS::BLAS IMPORTED UNKNOWN)
         set_target_properties(BLAS::BLAS PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${BLAS_LINKER_FLAGS}"
            IMPORTED_LOCATION "${BLAS_LIBRARIES}")
      endif()
   else()
      set(PAPILO_HAS_LUSOL 0)
   endif()
else()
   set(PAPILO_HAS_LUSOL 0)
   message(STATUS "LUSOL is disabled.")
endif()

# Determines whether -latomic must be added.
# TODO: Fix cmake/StdAtomic.cmake to create a target that can then be added to `libpapilo`.
include(cmake/StdAtomic.cmake)

# For `check_cxx_source_compiles` function.
include(CheckCXXSourceCompiles)
# For `check_cxx_compiler_flag` function.
include(CheckCXXCompilerFlag)

if(NOT MSVC)
   # Except for Windows we set the standard to C++14.
   check_cxx_compiler_flag(-std=c++14 HAVE_FLAG_STD_CXX14)
   if(HAVE_FLAG_STD_CXX14)
      set(CMAKE_REQUIRED_FLAGS "-std=c++14")
   else()
      check_cxx_compiler_flag(-std=c++1y HAVE_FLAG_STD_CXX1Y)
      if(HAVE_FLAG_STD_CXX1Y)
         set(CMAKE_REQUIRED_FLAGS "-std=c++1y")
      endif()
   endif()
endif()

# Includes for `check_cxx_source_compiles` call.
set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/src/papilo/external")
check_cxx_source_compiles(
   "#include \"ska/bytell_hash_map.hpp\"
   int main() { ska::bytell_hash_map<int,int> hashmap; (void)hashmap; return 0; }"
   PAPILO_HAS_BYTELL_HASHMAP)

if(PAPILO_BYTELL_HASHMAP_WORKS)
   set(PAPILO_USE_STANDARD_HASHMAP 0)
else()
   set(PAPILO_USE_STANDARD_HASHMAP 1)
endif()

##### TOOLS #####

set(CMAKE_CXX_STANDARD 14)
set(CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_EXTENSIONS OFF) # use -std=c++... instead of -std=gnu++...
# required for gcc if set(CMAKE_CXX_EXTENSIONS OFF) is used
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")

# if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
# endif()
if(MSVC)
   add_compile_options(/permissive)
endif()

#add_definitions ("-Wall")

# Disable fused floating point contraction to enhance reproducibility across compilers and
# architectures.
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
   add_compile_options(/fp:precise)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   add_compile_options(-fp-model=precise)
elseif((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
   add_compile_options(-ffp-contract=off)
endif()

##### SUBDIRECTORIES #####

# Main code
add_subdirectory(src)

# Create and install papilo binary unless we postpone it from within the scipoptsuite.
if(NOT PAPILO_NO_BINARIES)
   add_subdirectory(binaries)
endif()

# Add tests and include the subdirectory.
# TODO: Does this avoid PAPILO tests if included from scipoptsuite?!
if(BUILD_TESTING)
   enable_testing()
   add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif()

##### INSTALLATION #####

# Install license files.
install(FILES
   "${PROJECT_SOURCE_DIR}/LICENSE"
   DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/papilo"
   )
if(NOT PAPILO_HAS_SYSTEM_FMT)
   install(FILES
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/LICENSE.rst"
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/papilo/fmt"
      )
endif()
if(PAPILO_HAS_BYTELL_HASHMAP)
   install(FILES
      "${PROJECT_SOURCE_DIR}/src/papilo/external/ska/LICENSE.txt"
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/papilo/ska"
      )
endif()
if(PAPILO_HAS_LUSOL)
   install(FILES
      "${PROJECT_SOURCE_DIR}/src/papilo/external/lusol/LICENSE.md"
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/papilo/lusol"
      )
endif()

# Create papilo-config-version.cmake.
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_BINARY_DIR}/papilo-config-version.cmake"
   VERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}.${PAPILO_VERSION_PATCH}
   COMPATIBILITY SameMajorVersion
)

# Unless we build a package, we populate search directories for dependencies.
include(PackageConfigDependencies)

if(NOT BUILD_PACKAGE AND PAPILO_HAS_GMP AND GMP_DIR)
   package_config_dependency_dir(PAPILO_CONFIG_GMP_DIR GMP_FOUND GMP_DIR)
endif()
if(NOT BUILD_PACKAGE AND PAPILO_HAS_QUADMATH AND QUADMATH_DIR)
   package_config_dependency_dir(PAPILO_CONFIG_QUADMATH_DIR Quadmath_FOUND QUADMATH_DIR)
endif()
if(NOT BUILD_PACKAGE AND PAPILO_HAS_BOOST AND Boost_DIR)
   package_config_dependency_dir(PAPILO_CONFIG_BOOST_DIR Boost_FOUND Boost_DIR)
endif()
if(NOT BUILD_PACKAGE AND PAPILO_HAS_TBB AND TBB_DIR)
   package_config_dependency_dir(PAPILO_CONFIG_TBB_DIR TBB_FOUND TBB_DIR)
endif()

# Configure papilo-config.cmake for the build-tree:
# We add the CMake module path from the sources.
set(EXTRA_CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
configure_file(${PROJECT_SOURCE_DIR}/papilo-config.cmake.in
   "${CMAKE_BINARY_DIR}/papilo-config.cmake"
   @ONLY)

# Configure papilo-config.cmake for the installation.
# We add the current directory of the installed papilo-config.cmake.
set(EXTRA_CMAKE_MODULE_PATH "\${CMAKE_CURRENT_LIST_DIR}")
configure_file(${PROJECT_SOURCE_DIR}/papilo-config.cmake.in
   "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/papilo-config.cmake"
   @ONLY)

# Install CMake config and modules.
install(FILES
   "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/papilo-config.cmake"
   "${CMAKE_BINARY_DIR}/papilo-config-version.cmake"
   "$<$<BOOL:${PAPILO_HAS_GMP}>:${PROJECT_SOURCE_DIR}/cmake/FindGMP.cmake>"
   "$<$<BOOL:${PAPILO_HAS_QUADMATH}>:${PROJECT_SOURCE_DIR}/cmake/FindQuadmath.cmake>"
   "$<$<BOOL:${PAPILO_HAS_TBB}>:${PROJECT_SOURCE_DIR}/cmake/FindTBB.cmake>"
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/papilo"
   )

# Install the papilo-targets export set.
install(EXPORT papilo-targets
   NAMESPACE PAPILO::
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/papilo"
   )

# Export targets for the build-tree.
export(EXPORT papilo-targets
   NAMESPACE PAPILO::
   FILE "${CMAKE_BINARY_DIR}/papilo-targets.cmake")

# If we include this directory then we mark PaPILO as found.
set(PAPILO_FOUND 1)
