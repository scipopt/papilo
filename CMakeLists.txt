cmake_minimum_required (VERSION 3.11.0)

set(PAPILO_VERSION_MAJOR 2)
set(PAPILO_VERSION_MINOR 4)
set(PAPILO_VERSION_PATCH 1)
set(PAPILO_API_VERSION 3)

project(papilo VERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}.${PAPILO_VERSION_PATCH}.0)

set(BOOST_MIN_VERSION 1.65)
set(CMAKE_CXX_STANDARD 14)
set(CXX_STANDARD_REQUIRED ON)

# Each directory used during linking is appended to the rpath of installed libraries.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# set(CMAKE_CXX_EXTENSIONS OFF) # use -std=c++... instead of -std=gnu++...
# required for gcc if set(CMAKE_CXX_EXTENSIONS OFF) is used
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")

# if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
# endif()
if(MSVC)
   add_compile_options(/permissive)
endif()

# Sets installation directories according to GNU coding standards.
include(GNUInstallDirs)
#add_definitions ("-Wall")

# Location for storing CMake files.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules/)

# Functions without PAPILOLIB_EXPORT will not be accessible from outside.
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# Even static libs have position-independent code by default, and can be linked into shared libs.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Disable fused floating point contraction to enhance reproducibility across compilers and
# architectures.
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
   add_compile_options(/fp:precise)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   add_compile_options(-fp-model=precise)
elseif((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
   add_compile_options(-ffp-contract=off)
endif()

option(GMP "should gmp be linked" ON)
option(QUADMATH "should quadmath library be used" ON)
if(MSVC)
   option(LUSOL "should LUSOL package be enabled" OFF)
else()
   option(LUSOL "should LUSOL package be enabled" ON)
endif()

option(TBB "should TBB be linked if found" ON)
option(TBB_DOWNLOAD "should TBB be downloaded" OFF)
option(INSTALL_TBB "should the TBB library be installed" OFF)

option(SOPLEX "should SoPlex LP solver be linked if found" ON)
option(SCIP "should SCIP solver be linked if found" ON)
option(HIGHS "should HiGHS LP solver be linked" OFF)
option(GLOP "should ortools solver be linked" OFF)
option(GUROBI "should gurobi solver be linked" OFF)

option(BUILD_TESTING "should tests be enabled and built" ON)

# Make 'Release' the default build type.
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Search for Boost: headers are required and the other components are needed for the binaries.
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED OPTIONAL_COMPONENTS iostreams program_options serialization)
if(Boost_IOSTREAMS_FOUND)
   set(PAPILO_HAVE_BOOST_IOSTREAMS 1)
else()
   set(PAPILO_HAVE_BOOST_IOSTREAMS 0)
endif()
if(Boost_PROGRAM_OPTIONS_FOUND)
   set(PAPILO_COMMAND_LINE_AVAILABLE 1)
else()
   set(PAPILO_COMMAND_LINE_AVAILABLE 0)
endif()
if(Boost_SERIALIZATION_FOUND)
   set(PAPILO_SERIALIZATION_AVAILABLE 1)
endif()
if(NOT Boost_IOSTREAMS_FOUND OR NOT Boost_PROGRAM_OPTIONS_FOUND OR NOT Boost_SERIALIZATION_FOUND)
   message("WARNING: Boost iostreams, program_options or serialization not found. PaPILO can only be used as a library!")
endif()

# Search for GMP.
if(GMP)
   find_package(GMP)
endif()
if(GMP_FOUND)
   set(PAPILO_HAVE_GMP 1)
   set(GMP_VERSION "${GMP_VERSION}")
else()
   set(PAPILO_HAVE_GMP 0)
   set(GMP_INCLUDE_DIRS "")
   set(GMP_LIBRARIES "")
endif()

# Search for Quadmath.
if(QUADMATH)
   find_package(Quadmath)
endif()
if(Quadmath_FOUND)
   set(PAPILO_HAVE_FLOAT128 1)
else()
   set(PAPILO_HAVE_FLOAT128 0)
   set(Quadmath_INCLUDE_DIRS "")
endif()

# Search for Gurobi.
if(GUROBI)
   find_package(GUROBI REQUIRED)
endif()

# Search for GLOP.
if(GLOP)
   find_package(GLOP REQUIRED)
endif()

# Search for git.
find_program(GIT git)
if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
   # If git is available and this is a repository, we query the git hash.
   execute_process(
      COMMAND ${GIT} rev-parse --short HEAD
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE PAPILO_GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
   # Otherwise we extract the githash in case was defined in papilo/Config.hpp.
   file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/src/papilo/Config.hpp githash_define REGEX "define PAPILO_GITHASH .*")
   if(githash_define)
      string(REGEX MATCH "\\\"(.*)\\\"" _matched ${githash_define})
      if(_matched)
         set(PAPILO_GITHASH ${CMAKE_MATCH_1})
      endif()
   endif()
endif()

# Output the githash.
if(PAPILO_GITHASH)
   message(STATUS "Git hash: ${PAPILO_GITHASH}")
   set(PAPILO_GITHASH_AVAILABLE 1)
else()
   set(PAPILO_GITHASH_AVAILABLE 0)
endif()

# Target `papilo` is an interface with all the dependencies.
add_library(papilo INTERFACE)
if(ROUNDINGSAT_DIR)
   set(ROUNDINGSAT_INCLUDE_DIRS ${ROUNDINGSAT_DIR}/src/)
endif()
target_include_directories(papilo INTERFACE
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
   $<INSTALL_INTERFACE:include>
   ${GMP_INCLUDE_DIRS}
   ${Quadmath_INCLUDE_DIRS}
   ${ROUNDINGSAT_INCLUDE_DIRS})

# Add the Boost dependencies that were found.
if(PAPILO_HAVE_BOOST_IOSTREAMS)
   target_link_libraries(papilo INTERFACE Boost::iostreams)
endif()
if(PAPILO_COMMAND_LINE_AVAILABLE)
   target_link_libraries(papilo INTERFACE Boost::program_options)
endif()
if(PAPILO_SERIALIZATION_AVAILABLE)
   target_link_libraries(papilo INTERFACE Boost::serialization)
endif()
if(Quadmath_FOUND)
   target_link_libraries(papilo INTERFACE Quadmath::quadmath)
endif()

# Search for TBB.
if(TBB)
   set(PAPILO_TBB on)
   # On windows we do not build TBB.
   if(WIN32)
      find_package(TBB 2018 COMPONENTS tbb tbbmalloc REQUIRED)
   else()
      if(TBB_DOWNLOAD)
         # Utilities to automatically download and build missing dependencies.
         include(cmake/Dependencies.cmake)
         find_or_download_package(
            NAME TBB
            VERSION 2018
            URL https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.13.0.tar.gz
            COMPONENTS tbb tbbmalloc
            URL_HASH SHA256=3ad5dd08954b39d113dc5b3f8a8dc6dc1fd5250032b7c491eb07aed5c94133e1
            CONFIGURE_ARGS -DTBB_TEST=OFF -DTBB_EXAMPLES=OFF -DTBB4PY_BUILD=OFF)
      else()
         # Otherwise we try to find TBB.
         find_package(TBB 2018 COMPONENTS tbb tbbmalloc)
      endif()
      if(NOT TBB_FOUND)
         message(FATAL_ERROR "TBB was not found. Please specify its location (-DTBB_DIR=...), enable auto-download of TBB (-DTBB_DOWNLOAD=on) or disable TBB (-DTBB=off).")
      else()
         message(STATUS "Found TBB: ${TBB_DIR}") # TODO: Extract and print version information.
      endif()
   endif()

   # Add this to the `papilo` target.
   target_link_libraries(papilo INTERFACE TBB::tbb Threads::Threads $<$<PLATFORM_ID:Linux>:rt>)
else()
   set(PAPILO_TBB off)
   message(WARNING "Disabling TBB results in a dramatically degradation in the performance of PaPILO!")
endif()

# Search for Threads.
find_package(Threads REQUIRED)

# TODO: Fix FindGMP.cmake and FindGLOP to create targets instead of variables.
target_link_libraries(papilo
   INTERFACE ${GMP_LIBRARIES} ${GLOP_LIBRARIES})

# Determines whether -latomic must be added.
# TODO: Fix cmake/StdAtomic.cmake to create a target that can then be added to `papilo`.
include(cmake/StdAtomic.cmake)

# For `write_basic_package_version_file` function.
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
   ${CMAKE_BINARY_DIR}/papilo-config-version.cmake
   VERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}.${PAPILO_VERSION_PATCH}
   COMPATIBILITY SameMajorVersion
)

# Specific compilation flags.
if(MSVC)
   target_compile_definitions(papilo INTERFACE BOOST_ALL_NO_LIB)
   target_compile_options(papilo INTERFACE /bigobj)
   target_compile_options(papilo INTERFACE ${TBB_CXX_STD_FLAG})
else()
   #TODO: configure warning level https://stackoverflow.com/questions/2368811/how-to-set-warning-level-in-cmake
   target_compile_options(papilo INTERFACE ${TBB_CXX_STD_FLAG} -Wno-shadow -Wall)
endif()

# Search for BLAS and build clusol.
if(LUSOL)
   include(CheckLanguage)
   check_language(Fortran)
   find_package(BLAS)
endif()
if(CMAKE_Fortran_COMPILER AND BLAS_FOUND)
   enable_language(Fortran)

   add_library(clusol STATIC
      src/papilo/external/lusol/src/lusol6b.f
      src/papilo/external/lusol/src/lusol7b.f
      src/papilo/external/lusol/src/lusol8b.f
      src/papilo/external/lusol/src/lusol_util.f
      src/papilo/external/lusol/src/lusol.f90
      src/papilo/external/lusol/src/clusol.c
      src/papilo/external/lusol/src/lusol_precision.f90)
   set_property(TARGET clusol PROPERTY POSITION_INDEPENDENT_CODE ON)
   target_include_directories(clusol PRIVATE src/papilo/external/lusol)
   target_link_libraries(clusol ${BLAS_LIBRARIES})
   target_link_libraries(papilo INTERFACE clusol)
   set(PAPILO_HAVE_LUSOL 1)
   set(CLUSOL_TARGET clusol)
else()
   set(PAPILO_HAVE_LUSOL 0)
   set(CLUSOL_TARGET "")
   message(STATUS "LUSOL is not built")
endif()

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
set(CMAKE_REQUIRED_INCLUDES ${PROJECT_SOURCE_DIR}/src/papilo/external)
check_cxx_source_compiles(
   "#include \"ska/bytell_hash_map.hpp\"
   int main() { ska::bytell_hash_map<int,int> hashmap; (void)hashmap; return 0; }"
   PAPILO_BYTELL_HASHMAP_WORKS )

if(PAPILO_BYTELL_HASHMAP_WORKS)
   set(PAPILO_USE_STANDARD_HASHMAP 0)
else()
   set(PAPILO_USE_STANDARD_HASHMAP 1)
endif()

# Create target `papilo-core` which is the static library.
add_library(papilo-core STATIC
   src/papilo/core/VariableDomains.cpp
   src/papilo/core/SparseStorage.cpp
   src/papilo/core/ConstraintMatrix.cpp
   src/papilo/core/ProblemUpdate.cpp
   src/papilo/core/Presolve.cpp
   src/papilo/core/postsolve/PostsolveStorage.cpp
   src/papilo/core/postsolve/Postsolve.cpp
   src/papilo/core/ProbingView.cpp
   src/papilo/presolvers/CoefficientStrengthening.cpp
   src/papilo/presolvers/ConstraintPropagation.cpp
   src/papilo/presolvers/DominatedCols.cpp
   src/papilo/presolvers/DualFix.cpp
   src/papilo/presolvers/DualInfer.cpp
   src/papilo/presolvers/FixContinuous.cpp
   src/papilo/presolvers/FreeVarSubstitution.cpp
   src/papilo/presolvers/ImplIntDetection.cpp
   src/papilo/presolvers/ParallelColDetection.cpp
   src/papilo/presolvers/ParallelRowDetection.cpp
   src/papilo/presolvers/Probing.cpp
   src/papilo/presolvers/SimpleProbing.cpp
   src/papilo/presolvers/SimpleSubstitution.cpp
   src/papilo/presolvers/SingletonCols.cpp
   src/papilo/presolvers/SingletonStuffing.cpp
   src/papilo/presolvers/Sparsify.cpp
   src/papilo/presolvers/SimplifyInequalities.cpp
   ${roundingsat_source_files}
  )

#
# interface function for setting common library properties
#
# TODO: Simplify with a generator expression.
function(setLibProperties targetname outputname)
    set_target_properties(${targetname} PROPERTIES
        OUTPUT_NAME ${outputname}
        MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endfunction(setLibProperties)
if(MSVC)
    setLibProperties(papilo-core "libpapilo-core")
else()
   setLibProperties(papilo-core "papilo-core")
endif()


# Add dependency-target `papilo` to static library-target `papilo-core`.
target_link_libraries(papilo-core papilo)
target_compile_definitions(papilo-core PRIVATE PAPILO_USE_EXTERN_TEMPLATES)

# When including from scipoptsuite, this is set to postpone generation of binaries until SCIP is built.
if(NOT PAPILO_NO_BINARIES)
   add_subdirectory(binaries)
endif()

# Create build/papilo/CMakeConfig.hpp, which contains all definitions based on the configuration.
# TODO: Are quotes really needed?
configure_file("${PROJECT_SOURCE_DIR}/src/papilo/CMakeConfig.hpp.in"
   "${PROJECT_BINARY_DIR}/papilo/CMakeConfig.hpp")

# Add tests and include the subdirectory.
# TODO: Does this avoid PAPILO tests if included from scipoptsuite?!
if(BUILD_TESTING)
   enable_testing()
   add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif()

# Install the header files of PAPILO.
install(FILES
   ${PROJECT_BINARY_DIR}/papilo/CMakeConfig.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/Config.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo)

install(FILES
   ${CMAKE_BINARY_DIR}/papilo-config-version.cmake
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo)

# TODO: Any reason for ${PROJECT_SOURCE_DIR}/src as opposed to just src/ ?
install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Components.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ConstraintMatrix.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/MatrixBuffer.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Objective.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Presolve.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/PresolveMethod.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/PresolveOptions.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ProbingView.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Problem.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ProblemBuilder.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ProblemFlag.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ProblemUpdate.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Reductions.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/RowFlags.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/SingleRow.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Solution.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/SparseStorage.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Statistics.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/SymmetryStorage.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/VariableDomains.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/core)


install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/core/postsolve/BoundStorage.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/postsolve/PostsolveStorage.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/postsolve/Postsolve.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/postsolve/PostsolveStatus.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/postsolve/PostsolveType.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/postsolve/ReductionType.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/postsolve/SavedRow.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/core/postsolve)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/HighsInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/GlopInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/GurobiInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/ScipInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/SolverInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/SoplexInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/RoundingsatInterface.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/interfaces)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/io/BoundType.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/Message.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/MpsParser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/MpsWriter.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/OpbParser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/OpbWriter.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/ParseKey.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/Parser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/SolParser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/SolWriter.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/io)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Alloc.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Array.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/compress_vector.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/DependentRows.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Flags.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/fmt.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Hash.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/MultiPrecision.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Num.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/NumericalStatistics.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/PrimalDualSolValidation.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/OptionsParser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/VersionLogger.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/ParameterSet.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Signature.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/StableSum.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/String.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/tbb.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Timer.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Validation.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Vec.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/VectorUtils.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Wrappers.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/misc)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/verification/ArgumentType.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/verification/CertificateInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/verification/EmptyCertificate.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/verification/VeriPb.hpp
   DESTINATION include/papilo/verification)

# Add extended Euclidean code in case Boost version is below 1.70.
if((DEFINED Boost_VERSION) AND (Boost_VERSION VERSION_LESS 107000))
   install(FILES
      ${PROJECT_SOURCE_DIR}/src/papilo/misc/extended_euclidean.hpp
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/misc)
endif()

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/CoefficientStrengthening.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/ConstraintPropagation.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/DominatedCols.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/DualFix.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/DualInfer.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/FixContinuous.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/FreeVarSubstitution.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/ImplIntDetection.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/ParallelColDetection.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/ParallelRowDetection.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/Probing.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/SimpleProbing.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/SimpleSubstitution.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/SimplifyInequalities.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/SingletonCols.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/SingletonStuffing.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/presolvers/Sparsify.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/presolvers)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/chrono.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/color.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/compile.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/core.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/format-inl.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/format.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/locale.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/os.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/ostream.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/posix.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/printf.h
   ${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/ranges.h
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/external/fmt)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/external/pdqsort/pdqsort.h
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/external/pdqsort)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/external/ska/bytell_hash_map.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/external/ska/flat_hash_map.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/external/ska/unordered_map.hpp
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/external/ska)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/external/lusol/clusol.h
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/papilo/external/lusol)

# Configure the papilo-config.cmake for the build-tree:
# We add the CMake module path from the sources.
set(EXTRA_CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
configure_file(${PROJECT_SOURCE_DIR}/papilo-config.cmake.in
   "${CMAKE_BINARY_DIR}/papilo-config.cmake" @ONLY)

# Configure the papilo-config.cmake for the installation.
# We add the current directory of the installed papilo-config.cmake.
set(EXTRA_CMAKE_MODULE_PATH "\${CMAKE_CURRENT_LIST_DIR}")
configure_file(${PROJECT_SOURCE_DIR}/papilo-config.cmake.in
   "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/papilo-config.cmake" @ONLY)

# Install FindGMP.cmake module.
if(PAPILO_HAVE_GMP)
   install(FILES ${PROJECT_SOURCE_DIR}/cmake/Modules/FindGMP.cmake
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo)
endif()

# Install FindQuadmath.cmake module.
if(PAPILO_HAVE_FLOAT128)
   install(FILES ${PROJECT_SOURCE_DIR}/cmake/Modules/FindQuadmath.cmake
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo)
endif()

# Install the papilo-config.
install(FILES ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/papilo-config.cmake
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo)

# Install the FindTBB.cmake module.
install(FILES ${PROJECT_SOURCE_DIR}/cmake/Modules/FindTBB.cmake
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo)

# Install TBB if requested.
if(INSTALL_TBB AND CMAKE_SYSTEM_NAME STREQUAL "Windows")
   if(EXISTS "${TBB_tbb_DLL}")
      install(FILES "${TBB_tbb_DLL}" DESTINATION ${CMAKE_INSTALL_BINDIR})
      install(FILES "${TBB_tbb_DLL}" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo)
   endif()
   if(EXISTS "${TBB_tbbmalloc_DLL}")
      install(FILES "${TBB_tbbmalloc_DLL}" DESTINATION ${CMAKE_INSTALL_BINDIR})
   endif()
endif()

# Export targets for the build-tree.
export(TARGETS papilo papilo-core ${CLUSOL_TARGET} FILE "${CMAKE_BINARY_DIR}/papilo-targets.cmake")

# Install the dependencies.
if(EXISTS ${CMAKE_BINARY_DIR}/_deps)
   install(DIRECTORY ${CMAKE_BINARY_DIR}/_deps DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo/)
endif()

# Install targets.
install(TARGETS papilo papilo-core ${CLUSOL_TARGET} EXPORT papilo-targets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the papilo-targets.
install(EXPORT papilo-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/papilo)

