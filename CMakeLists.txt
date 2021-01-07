cmake_minimum_required (VERSION 3.5.0)
project(papilo VERSION 1.0.2.0)

set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_EXTENSIONS OFF)
# required for gcc if set(CMAKE_CXX_EXTENSIONS OFF) is used
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif()

# path to e.g. findGMP module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(GMP "should gmp be linked" ON)
option(QUADMATH "should quadmath library be used" ON)
if(MSVC)
   option(LUSOL "should LUSOL package be enabled" OFF)
else()
   option(LUSOL "should LUSOL package be enabled" ON)
endif()
option(HIGHS "should HiGHS LP solver be linked if found" ON)
option(SOPLEX "should SoPlex LP solver be linked if found" ON)
option(SCIP "should SCIP solver be linked if found" ON)

# make 'Release' the default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Boost REQUIRED)

if(GMP)
    find_package(GMP)
endif()
if(GMP_FOUND)
    set(PAPILO_HAVE_GMP 1)
else()
    set(GMP_INCLUDE_DIRS "")
    set(GMP_LIBRARIES "")
    set(PAPILO_HAVE_GMP 0)
endif()

if(QUADMATH)
   find_package(Quadmath)
endif()
if(Quadmath_FOUND)
   set(PAPILO_HAVE_FLOAT128 1)
   set(Quadmath_IMPORTED_TARGET Quadmath::quadmath)
else()
   set(PAPILO_HAVE_FLOAT128 0)
   set(Quadmath_INCLUDE_DIRS "")
   set(Quadmath_IMPORTED_TARGET "")
endif()

find_program(GIT git)

if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
   execute_process(
      COMMAND ${GIT} describe --always --dirty
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE PAPILO_GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REGEX REPLACE \"^.*-g\" \"\" PAPILO_GITHASH ${PAPILO_GITHASH})
else()
   file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/src/papilo/Config.hpp githash_define REGEX "define PAPILO_GITHASH .*")
   if(githash_define)
      string(REGEX MATCH "\\\"(.*)\\\"" _matched ${githash_define})
      if(_matched)
         set(PAPILO_GITHASH ${CMAKE_MATCH_1})
      endif()
   endif()
endif()

if(PAPILO_GITHASH)
   message(STATUS "Git hash: ${PAPILO_GITHASH}")
   set(PAPILO_GITHASH_AVAILABLE 1)
else()
   set(PAPILO_GITHASH_AVAILABLE 0)
endif()

add_library(papilo INTERFACE)
target_include_directories(papilo INTERFACE
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external>
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external/tbb/include>
   $<INSTALL_INTERFACE:include>
   $<INSTALL_INTERFACE:include/papilo-external>
   ${Boost_INCLUDE_DIRS}
   ${GMP_INCLUDE_DIRS}
   ${Quadmath_INCLUDE_DIRS})

find_package(TBB MODULE REQUIRED)
if(TBB_BUILT_STATIC_LIB)
   install(FILES ${TBB_STATIC_LIB} DESTINATION lib)
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
   file(TO_CMAKE_PATH "$ENV{TBBROOT}" TBB_ROOTDIR)
   if(CMAKE_CL_64)
      set(TBB_DLL "${TBB_ROOTDIR}/../redist/intel64/tbb/vc_mt/tbb.dll")
   else()
      set(TBB_DLL "${TBB_ROOTDIR}/../redist/ia32/tbb/vc_mt/tbb.dll")
   endif()
   if(EXISTS "${TBB_DLL}")
      install(FILES "${TBB_DLL}" DESTINATION bin)
   endif()
endif()

target_link_libraries(papilo INTERFACE ${TBB_LIBRARY} ${Quadmath_IMPORTED_TARGET} ${GMP_LIBRARIES})

if(MSVC)
   target_compile_definitions(papilo INTERFACE BOOST_ALL_NO_LIB)
   target_compile_options(papilo INTERFACE /bigobj)
endif()

if(LUSOL)
   include(CheckLanguage)
   check_language(Fortran)
   find_package(BLAS)
endif()
if(CMAKE_Fortran_COMPILER AND BLAS_FOUND)
   enable_language(Fortran)

   add_library(clusol STATIC
      external/lusol/src/lusol6b.f
      external/lusol/src/lusol7b.f
      external/lusol/src/lusol8b.f
      external/lusol/src/lusol_util.f
      external/lusol/src/lusol.f90
      external/lusol/src/clusol.c
      external/lusol/src/lusol_precision.f90)
   set_property(TARGET clusol PROPERTY POSITION_INDEPENDENT_CODE ON)
   target_include_directories(clusol PRIVATE external/lusol)
   target_link_libraries(clusol ${BLAS_LIBRARIES})
   target_link_libraries(papilo INTERFACE clusol)
   set(PAPILO_HAVE_LUSOL 1)
   set(CLUSOL_TARGET clusol)
else()
   set(PAPILO_HAVE_LUSOL 0)
   set(CLUSOL_TARGET "")
   message(STATUS "LUSOL is not built")
endif()

include(CheckCXXSourceCompiles)
include(CheckCXXCompilerFlag)

if(NOT MSVC)
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
set(CMAKE_REQUIRED_INCLUDES ${PROJECT_SOURCE_DIR}/external)
check_cxx_source_compiles(
   "#include \"ska/bytell_hash_map.hpp\"
   int main() { ska::bytell_hash_map<int,int> hashmap; (void)hashmap; return 0; }"
   PAPILO_BYTELL_HASHMAP_WORKS )

if( PAPILO_BYTELL_HASHMAP_WORKS )
   set(PAPILO_USE_STANDARD_HASHMAP 0)
else()
   set(PAPILO_USE_STANDARD_HASHMAP 1)
endif()

add_library(papilo-core STATIC
   src/papilo/core/VariableDomains.cpp
   src/papilo/core/SparseStorage.cpp
   src/papilo/core/ConstraintMatrix.cpp
   src/papilo/core/ProblemUpdate.cpp
   src/papilo/core/Presolve.cpp
   src/papilo/core/Postsolve.cpp
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
   src/papilo/presolvers/SimplifyInequalities.cpp)

target_link_libraries(papilo-core papilo)
target_compile_definitions(papilo-core PRIVATE PAPILO_USE_EXTERN_TEMPLATES)

if(NOT PAPILO_NO_BINARIES)
   add_subdirectory(binaries)
endif()

configure_file("${PROJECT_SOURCE_DIR}/src/papilo/CMakeConfig.hpp.in"
               "${PROJECT_BINARY_DIR}/papilo/CMakeConfig.hpp")

# add tests
enable_testing()
add_subdirectory(test)

# install the header files of papilo
install(FILES
   ${PROJECT_BINARY_DIR}/papilo/CMakeConfig.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/Config.hpp
   DESTINATION include/papilo)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Components.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ConstraintMatrix.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/MatrixBuffer.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Objective.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Postsolve.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Presolve.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/PresolveMethod.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/PresolveOptions.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ProbingView.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ProblemBuilder.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Problem.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/ProblemUpdate.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Reductions.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/RowFlags.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/SingleRow.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Solution.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/SparseStorage.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/Statistics.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/core/VariableDomains.hpp
   DESTINATION include/papilo/core)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/ScipInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/SolverInterface.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/interfaces/SoplexInterface.hpp
   DESTINATION include/papilo/interfaces)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/io/Message.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/MpsParser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/MpsWriter.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/SolParser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/io/SolWriter.hpp
   DESTINATION include/papilo/io)

install(FILES
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Alloc.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Array.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/compress_vector.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/DependentRows.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Flags.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/tbb.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/fmt.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Hash.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/KktChecker.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/KktCheckHelper.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/MultiPrecision.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Num.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/OptionsParser.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/ParameterSet.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Signature.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/StableSum.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/String.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Timer.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Vec.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/VectorUtils.hpp
   ${PROJECT_SOURCE_DIR}/src/papilo/misc/Wrappers.hpp
   DESTINATION include/papilo/misc)

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
   DESTINATION include/papilo/presolvers)

# install external library header files
install(DIRECTORY
   ${PROJECT_SOURCE_DIR}/external/tbb/include/tbb
   DESTINATION include/papilo-external)

install(FILES
   ${PROJECT_SOURCE_DIR}/external/fmt/chrono.h
   ${PROJECT_SOURCE_DIR}/external/fmt/color.h
   ${PROJECT_SOURCE_DIR}/external/fmt/compile.h
   ${PROJECT_SOURCE_DIR}/external/fmt/core.h
   ${PROJECT_SOURCE_DIR}/external/fmt/format.h
   ${PROJECT_SOURCE_DIR}/external/fmt/format-inl.h
   ${PROJECT_SOURCE_DIR}/external/fmt/locale.h
   ${PROJECT_SOURCE_DIR}/external/fmt/ostream.h
   ${PROJECT_SOURCE_DIR}/external/fmt/posix.h
   ${PROJECT_SOURCE_DIR}/external/fmt/printf.h
   ${PROJECT_SOURCE_DIR}/external/fmt/ranges.h
   ${PROJECT_SOURCE_DIR}/external/fmt/format.cc
   ${PROJECT_SOURCE_DIR}/external/fmt/posix.cc
   DESTINATION include/papilo-external/fmt)

install(FILES
   ${PROJECT_SOURCE_DIR}/external/pdqsort/pdqsort.h
   DESTINATION include/papilo-external/pdqsort)

install(FILES
   ${PROJECT_SOURCE_DIR}/external/ska/bytell_hash_map.hpp
   ${PROJECT_SOURCE_DIR}/external/ska/flat_hash_map.hpp
   ${PROJECT_SOURCE_DIR}/external/ska/unordered_map.hpp
   DESTINATION include/papilo-external/ska)

install(FILES
   ${PROJECT_SOURCE_DIR}/external/lusol/clusol.h
   DESTINATION include/papilo-external/lusol)

# configure the config file for the build tree
set(QUADMATH_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)
set(PAPILO_TBB_STATIC_LIB_LOCATION ${TBB_STATIC_LIB})
configure_file(${PROJECT_SOURCE_DIR}/papilo-config.cmake.in
  "${CMAKE_BINARY_DIR}/papilo-config.cmake" @ONLY)

# configure and install config file for installation
set(QUADMATH_MODULE_PATH "\${CMAKE_CURRENT_LIST_DIR}")
get_filename_component(TBB_STATIC_LIB_NAME "${TBB_STATIC_LIB}" NAME)
set(PAPILO_TBB_STATIC_LIB_LOCATION "\${CMAKE_CURRENT_LIST_DIR}/../../${TBB_STATIC_LIB_NAME}")
configure_file(${PROJECT_SOURCE_DIR}/papilo-config.cmake.in
    "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/papilo-config.cmake" @ONLY)

install(FILES ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/papilo-config.cmake
   DESTINATION lib/cmake/papilo)
if(PAPILO_HAVE_FLOAT128)
   install(FILES ${PROJECT_SOURCE_DIR}/cmake/Modules/FindQuadmath.cmake
      DESTINATION lib/cmake/papilo)
endif()

# export targets for build-tree linking
export(TARGETS papilo ${CLUSOL_TARGET} FILE "${CMAKE_BINARY_DIR}/papilo-targets.cmake")

# install targets
install(TARGETS papilo ${CLUSOL_TARGET} EXPORT papilo-targets ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(EXPORT papilo-targets DESTINATION lib/cmake/papilo)
