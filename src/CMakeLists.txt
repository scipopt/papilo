# Overall structure:
#
# Configure build-dependent header file.
# [Define common compilation settings.]
# Define targets along with properties, dependencies, includes, linked libraries, definitions.
# Define target aliases.
# [Add sanitizers.]
# Install headers.
# Install targets.
# Generate export header file.

# Create build/papilo/CMakeConfig.hpp, which contains all definitions based on the configuration.
if(BUILD_SHARED_LIBS)
   set(PAPILO_STATIC_LIB 0)
else()
   set(PAPILO_STATIC_LIB 1)
endif()
configure_file("${PROJECT_SOURCE_DIR}/src/papilo/CMakeConfig.hpp.in"
   "${PROJECT_BINARY_DIR}/papilo/CMakeConfig.hpp"
   @ONLY)

# Target `papilo` is an interface with all the dependencies.
add_library(libpapilo INTERFACE)

# Add include directories.
target_include_directories(libpapilo
   INTERFACE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>"
      "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/>"
      "$<INSTALL_INTERFACE:include/>")

# Add linked libraries.
target_link_libraries(libpapilo
   INTERFACE
      Boost::boost
      "$<$<BOOL:${PAPILO_HAS_GMP}>:GMP::GMP>"
      "$<$<BOOL:${PAPILO_HAS_QUADMATH}>:Quadmath::Quadmath>"
      "$<$<BOOL:${PAPILO_HAS_TBB}>:TBB::tbb>"
      "$<$<BOOL:${PAPILO_HAS_TBB}>:Threads::Threads>"
      "$<$<BOOL:${PAPILO_HAS_TBB}>:$<$<PLATFORM_ID:Linux>:rt>>"
      "$<$<BOOL:${PAPILO_HAS_SYSTEM_FMT}>:fmt::fmt>"
   )

# Add compilation flags.
target_compile_options(libpapilo
   INTERFACE
      $<$<BOOL:${MSVC}>:/bigobj>
      $<$<NOT:$<BOOL:${MSVC}>>:-Wno-shadow>
      $<$<NOT:$<BOOL:${MSVC}>>:-Wall>
   )

# Build clusol.
if(PAPILO_HAS_LUSOL)
   add_library(libclusol
      papilo/external/lusol/src/lusol6b.f
      papilo/external/lusol/src/lusol7b.f
      papilo/external/lusol/src/lusol8b.f
      papilo/external/lusol/src/lusol_util.f
      papilo/external/lusol/src/lusol.f90
      papilo/external/lusol/src/clusol.c
      papilo/external/lusol/src/lusol_precision.f90
      )
   set_target_properties(libclusol PROPERTIES
      PREFIX ""
      VERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}.${PAPILO_VERSION_PATCH}
      SOVERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}
      INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
      INSTALL_RPATH_USE_LINK_PATH TRUE
      )
   target_include_directories(libclusol PRIVATE "papilo/external/lusol")
   target_link_libraries(libclusol BLAS::BLAS)
   target_link_libraries(libpapilo INTERFACE libclusol)
endif()

# Create target `libpapilo-core` which is the main library.
if(PAPILO_HEADER_ONLY)
   add_library(libpapilo-core INTERFACE)
   # Add dependency-target `libpapilo` to library-target `libpapilo-core`.
   target_link_libraries(libpapilo-core INTERFACE libpapilo)
else()
   add_library(libpapilo-core
      papilo/core/VariableDomains.cpp
      papilo/core/SparseStorage.cpp
      papilo/core/ConstraintMatrix.cpp
      papilo/core/ProblemUpdate.cpp
      papilo/core/Presolve.cpp
      papilo/core/postsolve/PostsolveStorage.cpp
      papilo/core/postsolve/Postsolve.cpp
      papilo/core/ProbingView.cpp
      papilo/core/CliqueProbingView.cpp
      papilo/presolvers/CoefficientStrengthening.cpp
      papilo/presolvers/ConstraintPropagation.cpp
      papilo/presolvers/DominatedCols.cpp
      papilo/presolvers/DualFix.cpp
      papilo/presolvers/DualInfer.cpp
      papilo/presolvers/FixContinuous.cpp
      papilo/presolvers/GF2.cpp
      papilo/presolvers/FreeVarSubstitution.cpp
      papilo/presolvers/ImplIntDetection.cpp
      papilo/presolvers/ParallelColDetection.cpp
      papilo/presolvers/ParallelRowDetection.cpp
      papilo/presolvers/Probing.cpp
      papilo/presolvers/SimpleProbing.cpp
      papilo/presolvers/SimpleSubstitution.cpp
      papilo/presolvers/SingletonCols.cpp
      papilo/presolvers/SingletonStuffing.cpp
      papilo/presolvers/Sparsify.cpp
      papilo/presolvers/SimplifyInequalities.cpp
      papilo/presolvers/CliqueMerging.cpp
      ${roundingsat_source_files}
      )
   set_target_properties(libpapilo-core PROPERTIES
      PREFIX ""
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" # .dll
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" # .so
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" # .a
      VERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}.${PAPILO_VERSION_PATCH}
      SOVERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}
      INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
      INSTALL_RPATH_USE_LINK_PATH TRUE
      )
   # Add dependency-target `libpapilo` to library-target `libpapilo-core`.
   target_link_libraries(libpapilo-core
      PUBLIC
         libpapilo
      PRIVATE
         $<$<BOOL:{PAPILO_HAS_BOOST_SERIALIZATION}>:Boost::serialization>
      )
   # The PAPILO_BUILD_SHARED_LIB and PAPILO_BUILD_STATIC_LIB macros indicate that the shared / static
   # library is built.
   target_compile_definitions(libpapilo-core
      PRIVATE
         $<$<BOOL:${BUILD_SHARED_LIBS}>:PAPILO_BUILD_SHARED_LIB>
      )
endif()


# Deprecated papilo target. TODO: Remove in version 5.
add_library(papilo INTERFACE)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
   set_target_properties(papilo PROPERTIES
      DEPRECATION "papilo target is deprecated; use PAPILO::PAPILO instead.")
endif()
target_link_libraries(papilo INTERFACE libpapilo)
install(TARGETS papilo
   EXPORT papilo-targets
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Create aliases internally in case PaPILO is included via `add_subdirectory`.
add_library(PAPILO::PAPILO ALIAS libpapilo-core)
add_library(PAPILO::libpapilo-core ALIAS libpapilo-core)
add_library(PAPILO::libpapilo ALIAS libpapilo)
if(PAPILO_HAS_LUSOL)
   add_library(PAPILO::libclusol ALIAS libclusol)
endif()

# Install public header files.
foreach(header_file
      papilo/core/Components.hpp
      papilo/core/ConstraintMatrix.hpp
      papilo/core/MatrixBuffer.hpp
      papilo/core/Objective.hpp
      papilo/core/Presolve.hpp
      papilo/core/PresolveMethod.hpp
      papilo/core/PresolveOptions.hpp
      papilo/core/ProbingView.hpp
      papilo/core/CliqueProbingView.hpp
      papilo/core/Problem.hpp
      papilo/core/ProblemBuilder.hpp
      papilo/core/ProblemFlag.hpp
      papilo/core/ProblemUpdate.hpp
      papilo/core/Reductions.hpp
      papilo/core/RowFlags.hpp
      papilo/core/SingleRow.hpp
      papilo/core/Solution.hpp
      papilo/core/SparseStorage.hpp
      papilo/core/Statistics.hpp
      papilo/core/SymmetryStorage.hpp
      papilo/core/VariableDomains.hpp
      papilo/core/postsolve/BoundStorage.hpp
      papilo/core/postsolve/PostsolveStorage.hpp
      papilo/core/postsolve/Postsolve.hpp
      papilo/core/postsolve/PostsolveStatus.hpp
      papilo/core/postsolve/PostsolveType.hpp
      papilo/core/postsolve/ReductionType.hpp
      papilo/core/postsolve/SavedRow.hpp
      papilo/interfaces/HighsInterface.hpp
      papilo/interfaces/GlopInterface.hpp
      papilo/interfaces/GurobiInterface.hpp
      papilo/interfaces/ScipInterface.hpp
      papilo/interfaces/SolverInterface.hpp
      papilo/interfaces/SoplexInterface.hpp
      papilo/interfaces/RoundingsatInterface.hpp
      papilo/io/BoundType.hpp
      papilo/io/Message.hpp
      papilo/io/MpsParser.hpp
      papilo/io/MpsWriter.hpp
      papilo/io/OpbParser.hpp
      papilo/io/OpbWriter.hpp
      papilo/io/ParseKey.hpp
      papilo/io/Parser.hpp
      papilo/io/SolParser.hpp
      papilo/io/SolWriter.hpp
      papilo/misc/Alloc.hpp
      papilo/misc/Array.hpp
      papilo/misc/compress_vector.hpp
      papilo/misc/DependentRows.hpp
      papilo/misc/Flags.hpp
      papilo/misc/fmt.hpp
      papilo/misc/Hash.hpp
      papilo/misc/MultiPrecision.hpp
      papilo/misc/Num.hpp
      papilo/misc/NumericalStatistics.hpp
      papilo/misc/PrimalDualSolValidation.hpp
      papilo/misc/OptionsParser.hpp
      papilo/misc/VersionLogger.hpp
      papilo/misc/ParameterSet.hpp
      papilo/misc/Signature.hpp
      papilo/misc/StableSum.hpp
      papilo/misc/String.hpp
      papilo/misc/tbb.hpp
      papilo/misc/Timer.hpp
      papilo/misc/Validation.hpp
      papilo/misc/Vec.hpp
      papilo/misc/VectorUtils.hpp
      papilo/misc/Wrappers.hpp
      papilo/verification/ArgumentType.hpp
      papilo/verification/CertificateInterface.hpp
      papilo/verification/EmptyCertificate.hpp
      papilo/verification/VeriPb.hpp
      papilo/presolvers/CoefficientStrengthening.hpp
      papilo/presolvers/ConstraintPropagation.hpp
      papilo/presolvers/DominatedCols.hpp
      papilo/presolvers/DualFix.hpp
      papilo/presolvers/DualInfer.hpp
      papilo/presolvers/GF2.hpp
      papilo/presolvers/FixContinuous.hpp
      papilo/presolvers/FreeVarSubstitution.hpp
      papilo/presolvers/ImplIntDetection.hpp
      papilo/presolvers/ParallelColDetection.hpp
      papilo/presolvers/ParallelRowDetection.hpp
      papilo/presolvers/Probing.hpp
      papilo/presolvers/SimpleProbing.hpp
      papilo/presolvers/SimpleSubstitution.hpp
      papilo/presolvers/SimplifyInequalities.hpp
      papilo/presolvers/SingletonCols.hpp
      papilo/presolvers/SingletonStuffing.hpp
      papilo/presolvers/Sparsify.hpp
      papilo/presolvers/CliqueMerging.hpp
      papilo/external/fmt/chrono.h
      papilo/external/fmt/color.h
      papilo/external/fmt/compile.h
      papilo/external/fmt/core.h
      papilo/external/fmt/format-inl.h
      papilo/external/fmt/format.h
      papilo/external/fmt/locale.h
      papilo/external/fmt/os.h
      papilo/external/fmt/ostream.h
      papilo/external/fmt/posix.h
      papilo/external/fmt/printf.h
      papilo/external/fmt/ranges.h
      papilo/external/ska/bytell_hash_map.hpp
      papilo/external/ska/flat_hash_map.hpp
      papilo/external/ska/unordered_map.hpp
      papilo/Config.hpp
   )
   get_filename_component(directory "${header_file}" DIRECTORY)
   install(
      FILES "${PROJECT_SOURCE_DIR}/src/${header_file}"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${directory}")
endforeach()

# Install CMakeConfig.hpp
install(
   FILES "${PROJECT_BINARY_DIR}/papilo/CMakeConfig.hpp"
   DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/papilo/")

# Install the included FMT unless we use the system one.
if(NOT PAPILO_HAS_SYSTEM_FMT)
   install(FILES
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/chrono.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/color.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/compile.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/core.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/format-inl.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/format.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/locale.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/os.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/ostream.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/posix.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/printf.h"
      "${PROJECT_SOURCE_DIR}/src/papilo/external/fmt/ranges.h"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/papilo/external/fmt")
endif()

# Add extended Euclidean code in case Boost version is below 1.70.
if((DEFINED Boost_VERSION) AND (Boost_VERSION VERSION_LESS 107000))
   install(
      FILES "${PROJECT_SOURCE_DIR}/src/papilo/misc/extended_euclidean.hpp"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/papilo/misc/")
endif()

# Add LUSOL header only if supported.
if(PAPILO_HAS_LUSOL)
   install(
      FILES "${PROJECT_SOURCE_DIR}/src/papilo/external/lusol/clusol.h"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/papilo/external/lusol/")
endif()

# Install targets, and record them in papilo-targets export set.
install(TARGETS libpapilo libpapilo-core
   EXPORT papilo-targets
   RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
   ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
   LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
if(PAPILO_HAS_LUSOL)
   install(TARGETS libclusol
      EXPORT papilo-targets
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      )
endif()

