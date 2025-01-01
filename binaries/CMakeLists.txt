# Note: This file is included from papilo/CMakeLists.txt and scipoptsuite/CMakeLists.
#       All variables it relies on must be supplied by both these files.

include(GenerateExportHeader)

if(SCIP)
   find_package(SCIP)
endif()
if(SCIP_FOUND)
   message(STATUS "Found SCIP: ${SCIP_INCLUDE_DIRS}")
   set(PAPILO_HAVE_SCIP 1)
else()
   message(STATUS "Finding SCIP - not found")
   set(PAPILO_HAVE_SCIP 0)
endif()

if(GUROBI)
   find_package(GUROBI REQUIRED)
endif()
if(GUROBI_FOUND)
   message(STATUS "Found GUROBI: ${GUROBI_INCLUDE_DIRS}")
   set(PAPILO_HAVE_GUROBI 1)
   include_directories(${GUROBI_INCLUDE_DIRS})
else()
   message(STATUS "Finding GUROBI - not found")
   set(PAPILO_HAVE_GUROBI 0)
endif()

if(GLOP)
   find_package(GLOP REQUIRED)
endif()
if(GLOP_FOUND)
   message(STATUS "Found GLOP: ${GLOP_INCLUDE_DIRS}")
   set(PAPILO_HAVE_GLOP 1)
   include_directories(${GLOP_INCLUDE_DIRS})
else()
   message(STATUS "Finding GLOP - not found")
   set(PAPILO_HAVE_GLOP 0)
endif()


if(SOPLEX)
   find_package(SOPLEX)
endif()
if(SOPLEX_FOUND)
   message(STATUS "Found SOPLEX: ${SOPLEX_INCLUDE_DIRS}")
   set(PAPILO_HAVE_SOPLEX 1)
else()
   message(STATUS "Finding SOPLEX - not found")
   set(PAPILO_HAVE_SOPLEX 0)
endif()

if(HIGHS)
   find_package(HIGHS CONFIG REQUIRED)
endif()
if(HIGHS_FOUND)
   message(STATUS "Found HIGHS: ${HIGHS_INCLUDE_DIRS}")
   set(PAPILO_HAVE_HIGHS 1)
else()
   message(STATUS "Finding HIGHS - not found")
   set(PAPILO_HAVE_HIGHS 0)
endif()

if(ROUNDINGSAT_DIR)
   set(PAPILO_HAVE_ROUNDINGSAT 1)
endif()
set(CMAKE_CXX_STANDARD 14)
if(ROUNDINGSAT_DIR)
   message(STATUS "Found ROUNDINGSAT: ${ROUNDINGSAT_DIR}")
   set(CMAKE_CXX_STANDARD 17)
else()
   message(STATUS "Finding ROUNDINGSAT - not found")
   set(PAPILO_HAVE_ROUNDINGSAT 0)
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
find_package(Boost ${BOOST_MIN_VERSION} OPTIONAL_COMPONENTS iostreams program_options serialization)

include_directories(${SCIP_INCLUDE_DIRS} ${SOPLEX_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})

if(Boost_IOSTREAMS_FOUND)
   set(CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIRS})
   set(CMAKE_REQUIRED_LIBRARIES ${Boost_LIBRARIES})
   check_cxx_source_compiles(
      "#include <boost/iostreams/filter/gzip.hpp>
       int main() { auto decomp = boost::iostreams::gzip_decompressor(); (void)decomp; return 0; }"
      PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB )
   check_cxx_source_compiles(
      "#include <boost/iostreams/filter/bzip2.hpp>
       int main() { auto decomp = boost::iostreams::bzip2_decompressor(); (void)decomp; return 0; }"
      PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2 )
   unset(CMAKE_REQUIRED_INCLUDES)
   unset(CMAKE_REQUIRED_LIBRARIES)
endif()

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
         ${ROUNDINGSAT_DIR}/src/ConstrSimple.cpp
         ${ROUNDINGSAT_DIR}/src/LpSolver.cpp
         ${ROUNDINGSAT_DIR}/src/Solver.cpp
         ${ROUNDINGSAT_DIR}/src/parsing.cpp
         ${ROUNDINGSAT_DIR}/src/quit.cpp
         ${ROUNDINGSAT_DIR}/src/run.cpp
         ${ROUNDINGSAT_DIR}/src/SolverStructs.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/gpl_3_0.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/zib_academic.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/roundingsat.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/MIT.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/licenses.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/boost.cpp
         ${ROUNDINGSAT_DIR}/src/used_licenses/lgpl_3_0.cpp
   )
   endif()
   add_executable(papilo-executable ${CMAKE_CURRENT_LIST_DIR}/../src/papilo.cpp ${roundingsat_source_files})

   if(PAPILO_HAVE_ROUNDINGSAT)
      list(APPEND linked_solvers "Roundingsat")
      set(PAPILO_HAVE_ROUNDINGSAT 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAVE_ROUNDINGSAT)
   endif()
   if(PAPILO_HAVE_SCIP)
      list(APPEND linked_solvers "SCIP")
      set(PAPILO_HAVE_SCIP 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAVE_SCIP)
   endif()
   if(PAPILO_HAVE_GUROBI)
      list(APPEND linked_solvers "GUROBI")
      set(PAPILO_HAVE_GUROBI 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAVE_GUROBI)
   endif()
   if(PAPILO_HAVE_GLOP)
      list(APPEND linked_solvers "GLOP")
      set(PAPILO_HAVE_GLOP 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAVE_GLOP)
   endif()
   if(PAPILO_HAVE_SOPLEX)
      list(APPEND linked_solvers "SoPlex")
      set(PAPILO_HAVE_SOPLEX 1 PARENT_SCOPE)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAVE_SOPLEX)
   endif()
   if(PAPILO_HAVE_HIGHS)
      list(APPEND linked_solvers "HiGHS")
      set(PAPILO_HAVE_HIGHS 1 PARENT_SCOPE)
      find_package(OpenMP)
      target_compile_definitions(papilo-executable PRIVATE PAPILO_HAVE_HIGHS)
      target_link_libraries(papilo-executable highs::highs)
   endif()
   if(linked_solvers)
      message(STATUS "Linking solvers: ${linked_solvers}")
   else()
      message(STATUS "No solvers linked")
   endif()
   set_target_properties(papilo-executable PROPERTIES OUTPUT_NAME papilo RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
   target_link_libraries(papilo-executable papilo-core ${SOPLEX_LIBRARIES} ${SCIP_LIBRARIES} Boost::boost Boost::iostreams Boost::program_options Boost::serialization ${GUROBI_LIBRARIES})
   if(INSTALL_TBB)
      set_target_properties(papilo-executable PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/papilo/_deps/local/lib/)
   endif()
   target_compile_definitions(papilo-executable PRIVATE PAPILO_USE_EXTERN_TEMPLATES PAPILO_MPS_WRITER)
   install(TARGETS papilo-executable ARCHIVE DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

#   add_executable(convMPS EXCLUDE_FROM_ALL ${CMAKE_CURRENT_LIST_DIR}/../src/convMPS.cpp)
#   set_target_properties(convMPS PROPERTIES OUTPUT_NAME convMPS RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
#   target_link_libraries(convMPS papilo-core ${Boost_LIBRARIES})
else()
   message(WARNING "Executable of PaPILO is not built because Boost iostreams, serialization or program options is missing")
endif()

#TODO: consider case Gurobi and SCIP are build
if(PAPILO_HAVE_SCIP AND PAPILO_HAVE_SOPLEX)
   add_library(papilolib SHARED ${CMAKE_CURRENT_LIST_DIR}/../src/papilolib.cpp)
   set_target_properties(papilolib PROPERTIES OUTPUT_NAME papilo LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
   target_include_directories(papilolib INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
   target_link_libraries(papilolib papilo ${SOPLEX_PIC_LIBRARIES} ${SCIP_LIBRARIES})
   if( Boost_IOSTREAMS_FOUND )
      target_compile_definitions(papilolib PRIVATE PAPILO_MPS_WRITER)
      target_link_libraries(papilolib Boost::boost Boost::iostreams)
   endif()
   generate_export_header(papilolib)
endif()
