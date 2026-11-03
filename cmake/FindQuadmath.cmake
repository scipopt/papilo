# Finds the Quadmath C library or detects that it is incorporated in the compiler.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - Quadmath::Quadmath (interfaces Quadmath::libquadmath)
# - Quadmath::libquadmath (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - Quadmath_FOUND (bool indicating presence of Quadmath::libquadmath)
#
# == Hints ==
#
# QUADMATH_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)

# Search for the C library.
function(quadmath_find_main_lib result search_dir)
   message(DEBUG "    quadmath_find_main_lib(\"${search_dir}\") with quadmath_no_default_path=<${quadmath_no_default_path}>")
   find_library(quadmath_main_lib
      NAMES "quadmath"
      PATHS "${search_dir}" "${search_dir}/../lib"
      ${quadmath_no_default_path}
      )
   message(DEBUG "    -> ${quadmath_main_lib}, stored in ${result}")
   set(${result} "${quadmath_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(quadmath_validate result include_dir)
   message(DEBUG "  quadmath_validate(\"${the_include_dir}\") with quadmath_no_default_path=<${quadmath_no_default_path}>")

   # Check for library.
   quadmath_find_main_lib(quadmath_main_lib "${include_dir}")
   message(DEBUG "  -> quadmath_main_lib=${quadmath_main_lib}")
   if(NOT quadmath_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${quadmath_main_lib}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If Quadmath was already found before.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND TARGET Quadmath::Quadmath
      AND TARGET Quadmath::libquadmath
   )
   message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME}.")
   return()
endif()

# Main search function whose variables are local.
function(quadmath_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(quadmath_no_default_path NO_DEFAULT_PATH)
      else()
         set(quadmath_no_default_path)
      endif()

      find_path(quadmath_main_inc
         NAMES "quadmath.h"
         HINTS "${QUADMATH_DIR}" "$ENV{QUADMATH_DIR}"
         PATH_SUFFIXES "include"
         VALIDATOR quadmath_validate
         )

      message(DEBUG "With quadmath_no_default_path=<${quadmath_no_default_path}>, find_path returned quadmath_main_inc <${quadmath_main_inc}>.")
      if(quadmath_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${quadmath_main_inc}>.")
            quadmath_validate(quadmath_validation_result "${quadmath_main_inc}")
            if(quadmath_validation_result)
               break()
            else()
               unset(quadmath_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(quadmath_main_inc)
      quadmath_find_main_lib(quadmath_main_lib "${quadmath_main_inc}")
      set(quadmath_list_libs "${quadmath_main_lib}")
      if(quadmath_main_lib)

         # C library target
         if(NOT TARGET Quadmath::libquadmath)
            add_library(Quadmath::libquadmath UNKNOWN IMPORTED)
            set_target_properties(Quadmath::libquadmath PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${quadmath_main_inc}"
               IMPORTED_LOCATION "${quadmath_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET Quadmath::Quadmath)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(Quadmath::Quadmath INTERFACE IMPORTED)
            target_link_libraries(Quadmath::Quadmath INTERFACE Quadmath::libquadmath)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME}: ${quadmath_main_inc} and ${quadmath_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(quadmath_fail_reason ": unexpected failure (quadmath_main_lib=${quadmath_main_lib}, quadmath_cxx_inc=${quadmath_cxx_inc}, quadmath_cxx_lib=${quadmath_cxx_lib})")
      endif()
   else()
      # As a backup, test whether it just works; in some cases, header and library are in search paths of the compiler, but not in those of CMake.
      include(CMakePushCheckState)
      include(CheckCXXSourceCompiles)
   
      cmake_push_check_state(RESET)
      list(APPEND CMAKE_REQUIRED_LIBRARIES "quadmath")
      set(CMAKE_REQUIRED_QUIET TRUE)
      message(DEBUG "Running compilation test.")
      CHECK_CXX_SOURCE_COMPILES("
         #include <quadmath.h>
   
         int main(void){
            __float128 foo = ::sqrtq(123.456);
            (void) foo;
         }"
         use_directly)
   
      cmake_pop_check_state()
      message(DEBUG "Compilation test result: ${use_directly}")
      if (use_directly)
         add_library(Quadmath::libquadmath INTERFACE IMPORTED)
         set_target_properties(Quadmath::libquadmath PROPERTIES
            INTERFACE_LINK_LIBRARIES "quadmath")
         add_library(Quadmath::Quadmath INTERFACE IMPORTED)
         target_link_libraries(Quadmath::Quadmath INTERFACE Quadmath::libquadmath)

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME}: compiler integration.")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      endif()
 
      set(quadmath_fail_reason "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency Quadmath${quadmath_fail_reason}. Try specifying suitable QUADMATH_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find Quadmath${quadmath_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency Quadmath${quadmath_fail_reason}. Try specifying suitable QUADMATH_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
quadmath_find()

