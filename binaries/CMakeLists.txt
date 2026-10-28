# We include this directory either from papilo or separately from scipoptsuite.
# In the latter case, the PAPILO_HAS_* variables are local to papilo and thus
# not visible here. Hence, we search for PaPILO.
find_package(PAPILO CONFIG QUIET HINTS "${CMAKE_BINARY_DIR}")

# Search for SCIP.
set(PAPILO_HAS_SCIP 0)
if(SCIP)
   if(SCIP STREQUAL "AUTO")
      find_package(SCIP CONFIG HINTS "${SCIP_DIR}")
   else()
      find_package(SCIP REQUIRED CONFIG HINTS "${SCIP_DIR}")
   endif()
   if(SCIP_FOUND)
      set(PAPILO_HAS_SCIP 1)
   endif()
endif()

# Deal with legacy SCIP.
if(PAPILO_HAS_SCIP AND NOT TARGET SCIP::libscip)
   message(STATUS "Found SCIP version ${SCIP_VERSION}: ${SCIP_INCLUDE_DIRS} and ${SCIP_LIBRARIES}")
   if(TARGET libscip)
      message(WARNING "Defining target alias SCIP::libscip -> libscip. Consider using SCIP version at least 11.0!")
      add_library(SCIP::libscip ALIAS libscip)
   else()
      message(FATAL_ERROR "Expected CMake targets SCIP::libscip or libscip!")
   endif()
endif()

# Search for Gurobi.
if(GUROBI)
   if(GUROBI STREQUAL "AUTO")
      find_package(Gurobi COMPONENTS cxx)
   else()
      find_package(Gurobi REQUIRED COMPONENTS cxx)
   endif()
   if(Gurobi_FOUND)
      set(PAPILO_HAS_GUROBI 1)
   else()
      set(PAPILO_HAS_GUROBI 0)
   endif()
else()
   message(STATUS "Solver Gurobi disabled.")
   set(PAPILO_HAS_GUROBI 0)
endif()

# Search for GLOP.
set(PAPILO_HAS_GLOP 0)
if(GLOP)
   if(ortools_FOUND)
      message(STATUS "Already found suitable OR-Tools version ${ortools_VERSION}.")
   else()
      message(STATUS "Searching for OR-Tools.")
      # Indent messages for dependency search.
      list(APPEND CMAKE_MESSAGE_INDENT "  ")
      set(OLD_CMAKE_FIND_PACKAGE_PREFER_CONFIG ${CMAKE_FIND_PACKAGE_PREFER_CONFIG})
      set(CMAKE_FIND_PACKAGE_PREFER_CONFIG 1)
      if(GLOP STREQUAL "AUTO")
         find_package(ortools CONFIG HINTS "${ORTOOLS_DIR}")
      else()
         find_package(ortools REQUIRED CONFIG HINTS "${ORTOOLS_DIR}")
      endif()
      set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ${OLD_CMAKE_FIND_PACKAGE_PREFER_CONFIG})

      # Unindent messages again; equivalent to list(POP_BACK, CMAKE_MESSAGE_INDENT) which requires CMake 3.15.
      list(REVERSE CMAKE_MESSAGE_INDENT)
      list(REMOVE_AT CMAKE_MESSAGE_INDENT 0)
      list(REVERSE CMAKE_MESSAGE_INDENT)

      if(ortools_FOUND)
         message(STATUS "Found OR-Tools version ${ortools_VERSION}: ${ortools_CONFIG}")
      endif()
   endif()
   if(ortools_FOUND)
      set(PAPILO_HAS_GLOP 1)
   endif()
else()
   message(STATUS "GLOP is disabled.")
endif()

# Search for SoPlex.
if(SOPLEX)
   if(SOPLEX STREQUAL "AUTO")
      find_package(SOPLEX CONFIG HINTS "${SOPLEX_DIR}")
   else()
      find_package(SOPLEX REQUIRED CONFIG HINTS "${SOPLEX_DIR}")
   endif()
   if(SOPLEX_FOUND)
      # Deal with legacy SoPlex.
      if(NOT TARGET SOPLEX::SOPLEX)
         # In older versions, soplex-config.cmake does not tell whether SOPLEX was compiled with Boost or GMP.
         message(WARNING "Creating alias targets SOPLEX::SOPLEX for libsoplexshared; consider upgrading SoPlex.")
         add_library(SOPLEX::SOPLEX ALIAS libsoplexshared)
         set(CMAKE_REQUIRED_LIBRARIES libsoplexshared)
         if(NOT DEFINED SOPLEX_HAS_BOOST)
            check_cxx_source_compiles("
               #include <soplex/config.h>
               int main()
               {
               #if SOPLEX_HAS_BOOST || defined(SOPLEX_WITH_BOOST)
                  return 0;
               #else
                  error
               #endif
               }"
               SOPLEX_HAS_BOOST)
         endif()
         if(NOT DEFINED SOPLEX_HAS_GMP)
            check_cxx_source_compiles("
               #include <soplex/config.h>
               int main()
               {
               #if SOPLEX_HAS_GMP || defined(SOPLEX_WITH_GMP)
                  return 0;
               #else
                  error
               #endif
               }"
               SOPLEX_HAS_GMP)
         endif()
         unset(CMAKE_REQUIRED_LIBRARIES)
      endif()
      unset(_reason)
      if(PAPILO_HAS_GMP AND NOT SOPLEX_HAS_GMP)
         set(_reason "missing support for GMP")
      elseif(PAPILO_HAS_BOOST AND NOT SOPLEX_HAS_BOOST)
         set(_reason "missing support for Boost")
      elseif(NOT PAPILO_HAS_GMP AND SOPLEX_HAS_GMP)
         set(_reason "support for GMP which is absent for PaPILO")
      elseif(NOT PAPILO_HAS_BOOST AND SOPLEX_HAS_BOOST)
         set(_reason "support for Boost which is absent for PaPILO")
      endif()
      if(_reason)
         message(STATUS "Found SoPlex version ${SOPLEX_VERSION}, but disabled due to ${_reason}.")
         set(PAPILO_HAS_SOPLEX 0)
      else()
         if(SOPLEX_VERSION VERSION_LESS 9)
            message(STATUS "Found SoPlex version ${SOPLEX_VERSION}: ${SOPLEX_CONFIG}")
         endif()
         set(PAPILO_HAS_SOPLEX 1)
      endif()
   else()
      message(STATUS "Could NOT find SoPlex. Try specifying SOPLEX_DIR!")
      set(PAPILO_HAS_SOPLEX 0)
   endif()
else()
   message(STATUS "SoPlex LP solver is disabled.")
   set(PAPILO_HAS_SOPLEX 0)
endif()

# Search for HiGHS.
if(HIGHS)
   if(HIGHS STREQUAL "AUTO")
      find_package(HIGHS CONFIG HINTS "${HIGHS_DIR}")
   else()
      find_package(HIGHS REQUIRED CONFIG HINTS "${HIGHS_DIR}")
   endif()
   if(HIGHS_FOUND)
      message(STATUS "Found HiGHS version ${HIGHS_VERSION}: ${HIGHS_CONFIG}")
      set(PAPILO_HAS_HIGHS 1)
   else()
      set(PAPILO_HAS_HIGHS 0)
   endif()
else()
   set(PAPILO_HAS_HIGHS 0)
endif()

if(ROUNDINGSAT_DIR)
   set(PAPILO_HAS_ROUNDINGSAT 1)
endif()
set(CMAKE_CXX_STANDARD 14)
if(ROUNDINGSAT_DIR)
   message(STATUS "Found ROUNDINGSAT: ${ROUNDINGSAT_DIR}")
   set(CMAKE_CXX_STANDARD 20)
else()
   message(STATUS "Searching for ROUNDINGSAT - not found")
   set(PAPILO_HAS_ROUNDINGSAT 0)
endif()

# Enable the Fortran language support when the clusol target is included.
# Even though this was done in the papilo CMake to create the target, it must be
# done here as well, since the language support is restricted to the calling directory
# and subdirectories, and since the papilo/binaries CMake can be called from a
# CMake which doesn't have Fortran language enabled, Fortran might not be enabled here yet.
if(TARGET clusol)
   include(CheckLanguage)
   check_language(Fortran)
   enable_language(Fortran)
endif()

# Search again since we cannot rely on papilo/CMakeLists.txt having searched for Boost already.
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED OPTIONAL_COMPONENTS iostreams program_options serialization)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

if(Boost_IOSTREAMS_FOUND AND Boost_SERIALIZATION_FOUND AND Boost_PROGRAM_OPTIONS_FOUND)
#   add_executable(duplicates  ${CMAKE_CURRENT_LIST_DIR}/../src/duplicates.cpp)
#   set_target_properties(duplicates PROPERTIES OUTPUT_NAME duplicates RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
#   target_link_libraries(duplicates papilo-core ${Boost_LIBRARIES})
#   target_compile_definitions(duplicates PRIVATE PAPILO_USE_EXTERN_TEMPLATES PAPILO_MPS_WRITER)

   if(ROUNDINGSAT_DIR)
      set(
         roundingsat_source_files
         ${ROUNDINGSAT_DIR}/src/Constr.cpp
         ${ROUNDINGSAT_DIR}/src/ConstrExp.cpp
         ${ROUNDINGSAT_DIR}/src/ConstrExpPool.cpp
         ${ROUNDINGSAT_DIR}/src/ConstrSimple.cpp
         ${ROUNDINGSAT_DIR}/src/Env.cpp
         ${ROUNDINGSAT_DIR}/src/Logger.cpp
         ${ROUNDINGSAT_DIR}/src/LpSolver.cpp
         ${ROUNDINGSAT_DIR}/src/ProofBuffer.cpp
         ${ROUNDINGSAT_DIR}/src/Solver.cpp
         ${ROUNDINGSAT_DIR}/src/SolverStructs.cpp
         ${ROUNDINGSAT_DIR}/src/Stats.cpp
         ${ROUNDINGSAT_DIR}/src/objrewrite.cpp
         ${ROUNDINGSAT_DIR}/src/optimize.cpp
         ${ROUNDINGSAT_DIR}/src/parsing.cpp
         ${ROUNDINGSAT_DIR}/src/quit.cpp
         ${ROUNDINGSAT_DIR}/src/run.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/gpl_3_0.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/apache_2.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/roundingsat.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/MIT.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/licenses.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/boost.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/lgpl_3_0.cpp
   )
   endif()
   add_executable(papilo-executable ${CMAKE_CURRENT_LIST_DIR}/../src/papilo.cpp ${roundingsat_source_files})

   if(PAPILO_HAS_ROUNDINGSAT)
      list(APPEND linked_solvers "Roundingsat")
      set(PAPILO_HAS_ROUNDINGSAT 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAS_ROUNDINGSAT=1)
      target_include_directories(papilo-executable PRIVATE ${ROUNDINGSAT_DIR}/src)
   endif()
   if(PAPILO_HAS_SCIP)
      list(APPEND linked_solvers "SCIP")
      set(PAPILO_HAS_SCIP 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAS_SCIP=1)
      target_link_libraries(papilo-executable SCIP::libscip)
   endif()
   if(PAPILO_HAS_GUROBI)
      list(APPEND linked_solvers "Gurobi")
      set(PAPILO_HAS_GUROBI 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAS_GUROBI=1)
      target_link_libraries(papilo-executable Gurobi::libgurobixx)
   endif()
   if(PAPILO_HAS_GLOP)
      list(APPEND linked_solvers "GLOP")
      set(PAPILO_HAS_GLOP 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAS_GLOP=1)
      target_link_libraries(papilo-executable ortools::ortools)
   endif()
   if(PAPILO_HAS_SOPLEX)
      list(APPEND linked_solvers "SoPlex")
      set(PAPILO_HAS_SOPLEX 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAS_SOPLEX=1)
      target_link_libraries(papilo-executable SOPLEX::SOPLEX)
   endif()
   if(PAPILO_HAS_HIGHS)
      list(APPEND linked_solvers "HiGHS")
      set(PAPILO_HAS_HIGHS 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAS_HIGHS=1)
      target_link_libraries(papilo-executable highs::highs)
   endif()
   if(linked_solvers)
      message(STATUS "Linking solvers: ${linked_solvers}")
   else()
      message(STATUS "No solvers linked")
   endif()
   set_target_properties(papilo-executable PROPERTIES
      OUTPUT_NAME papilo
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
      INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
      INSTALL_RPATH_USE_LINK_PATH TRUE
      )
   target_link_libraries(papilo-executable libpapilo-core Boost::boost Boost::iostreams Boost::program_options Boost::serialization)
   target_compile_definitions(papilo-executable PRIVATE PAPILO_MPS_WRITER)
   if(MSVC)
       install(FILES $<TARGET_PDB_FILE:papilo-executable> DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)
   endif()

   # Add papilo-executable to the papilo-targets export set.
   install(TARGETS papilo-executable
      EXPORT papilo-targets
      ARCHIVE
         DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY
         DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME
         DESTINATION ${CMAKE_INSTALL_BINDIR})
else()
   message(WARNING "Executable of PaPILO is not built because Boost iostreams, serialization or program options is missing")
endif()

include(GenerateExportHeader)
if(PAPILO_HAS_SCIP AND PAPILO_HAS_SOPLEX)
   if(BUILD_SHARED_LIBS AND (DEFINED SCIP_HAS_PIC) AND (NOT SCIP_HAS_PIC))
      message(WARNING "Disabling shared C library since SCIP was built without position-independent code!")
   elseif(BUILD_SHARED_LIBS AND (DEFINED SOPLEX_HAS_PIC) AND (NOT SOPLEX_HAS_PIC))
      message(WARNING "Disabling shared C library since SoPlex was built without position-independent code!")
   else()
      add_library(libpapilo_c "${CMAKE_CURRENT_LIST_DIR}/../src/papilolib.cpp")
      add_library(PAPILO::libpapilo_c ALIAS libpapilo_c)
      set_target_properties(libpapilo_c PROPERTIES
         OUTPUT_NAME papilo_c
         LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
         VERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}.${PAPILO_VERSION_PATCH}
         SOVERSION ${PAPILO_VERSION_MAJOR}.${PAPILO_VERSION_MINOR}
         INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
         INSTALL_RPATH_USE_LINK_PATH TRUE
         )
      target_include_directories(libpapilo_c
         INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
         )
      target_link_libraries(libpapilo_c
         libpapilo SCIP::libscip
         SOPLEX::SOPLEX
         $<$<BOOL:${PAPILO_HAS_BOOST_IOSTREAMS}>:Boost::boost>
         $<$<BOOL:${PAPILO_HAS_BOOST_IOSTREAMS}>:Boost::iostreams>
         )
      target_compile_definitions(libpapilo_c
         PRIVATE
            $<$<BOOL:${PAPILO_HAS_BOOST_IOSTREAMS}>:PAPILO_MPS_WRITER>
         )
      generate_export_header(libpapilo_c
         EXPORT_FILE_NAME "papilolib_export.h"
         EXPORT_MACRO_NAME "PAPILOLIB_EXPORT"
         )

      # Add libpapilo_c to the papilo-targets export set.
      install(TARGETS libpapilo_c
         EXPORT papilo-targets
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         )
   endif()
endif()
