include(FindPackageHandleStandardArgs)

# Define search paths based on user input and environment variables
set(TBB_SEARCH_DIR ${TBB_LIBRARY_DIR} ${TBB_ROOT_DIR} $ENV{TBB_INSTALL_DIR} $ENV{TBBROOT})

# for windows add additional default search paths
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
   set(TBB_DEFAULT_SEARCH_DIR  "C:/Program Files/Intel/TBB"
                               "C:/Program Files (x86)/Intel/TBB")

   if (CMAKE_CL_64)
      list(APPEND TBB_LIB_PATH_SUFFIXES lib/intel64/vc14)
      list(APPEND TBB_LIB_PATH_SUFFIXES bin/intel64/vc14)
   else ()
      list(APPEND TBB_LIB_PATH_SUFFIXES lib/ia32/vc14)
      list(APPEND TBB_LIB_PATH_SUFFIXES bin/ia32/vc14)
   endif ()

endif()

# try to find the tbb library in the system
find_library(TBB_LIBRARY
            NAMES tbb libtbb.so.2
            HINTS ${TBB_SEARCH_DIR}
            PATHS ${TBB_DEFAULT_SEARCH_DIR}
            PATH_SUFFIXES ${TBB_LIB_PATH_SUFFIXES})
set(TBB_BUILT_STATIC_LIB 0)

# if the library was not found try to build a static library from source
if(NOT TBB_LIBRARY)
   include(${CMAKE_CURRENT_LIST_DIR}/../../external/tbb/cmake/TBBBuild.cmake)
   tbb_build(TBB_ROOT ${CMAKE_CURRENT_LIST_DIR}/../../external/tbb CONFIG_DIR TBB_DIR MAKE_ARGS extra_inc=big_iron.inc tbb_build_dir=${CMAKE_CURRENT_BINARY_DIR} tbb_build_prefix=tbb)
   if(TBB_DIR)
      # building was successful
      find_library(TBB_LIBRARY
         NAMES tbb
         HINTS
            ${CMAKE_CURRENT_BINARY_DIR}/tbb_release
            ${CMAKE_CURRENT_BINARY_DIR}/tbb_debug)
      if(TBB_LIBRARY)
         # create an imported target that also links the thread libraries via its interface
         set(TBB_BUILT_STATIC_LIB 1)
         find_package( Threads )
         list(APPEND linkinterface "${CMAKE_THREAD_LIBS_INIT}")
         find_library(LIBRT rt)
         if(LIBRT)
            list(APPEND linkinterface "${LIBRT}")
         endif()
         add_library(tbb STATIC IMPORTED)
         set_target_properties(tbb PROPERTIES
            IMPORTED_LOCATION "${TBB_LIBRARY}"
            INTERFACE_LINK_LIBRARIES "${linkinterface}")
         # add variable to hold the static lib so that it can be installed
         set(TBB_STATIC_LIB ${TBB_LIBRARY})
         set(TBB_LIBRARY tbb)
      endif()
   endif()
endif()

find_package_handle_standard_args(TBB REQUIRED_VARS TBB_LIBRARY)
