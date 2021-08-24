include(FindPackageHandleStandardArgs)

# Define search paths based on user input and environment variables
set(TBB_SEARCH_DIR ${TBB_LIBRARY_DIR} ${TBB_ROOT_DIR} ${TBB_DIR} $ENV{TBB_INSTALL_DIR} $ENV{TBBROOT})

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

if(APPLE)
    foreach (i tbb tbb@2020 tbb@2020_U1 tbb@2020_U2 tbb@2020_U3 tbb@2020_U3_1)
        foreach (j 2020_U1 2020_U2 2020_U3 2020_U3_1)
            list(APPEND TBB_LIB_PATH_SUFFIXES ${i}/${j}/lib)
        endforeach()
    endforeach()
    find_library(TBB_LIBRARY
                NAMES libtbb.dylib
                HINTS ${TBB_SEARCH_DIR}
                PATHS /usr/local/Cellar/
                NO_DEFAULT_PATH
                PATH_SUFFIXES ${TBB_LIB_PATH_SUFFIXES})
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    find_library(TBB_LIBRARY
                NAMES tbb
                HINTS ${TBB_SEARCH_DIR}
                PATHS ${TBB_DEFAULT_SEARCH_DIR}
                PATH_SUFFIXES ${TBB_LIB_PATH_SUFFIXES})
else()
    find_library(TBB_LIBRARY
                NAMES libtbb.so.2
                HINTS ${TBB_SEARCH_DIR}
                PATHS ${TBB_DEFAULT_SEARCH_DIR}
                PATH_SUFFIXES ${TBB_LIB_PATH_SUFFIXES})
endif()

set(TBB_BUILT_STATIC_LIB 0)

# if the library was not found try to build a static library from source
if((NOT TBB_LIBRARY) AND (NOT APPLE))
   include(${CMAKE_CURRENT_LIST_DIR}/../../external/tbb/cmake/TBBBuild.cmake)
   tbb_build(TBB_ROOT ${CMAKE_CURRENT_LIST_DIR}/../../external/tbb CONFIG_DIR TBB_DIR MAKE_ARGS extra_inc=big_iron.inc tbb_build_dir=${CMAKE_CURRENT_BINARY_DIR} tbb_build_prefix=tbb)
   if(TBB_DIR)
      # building was successful
      find_library(TBB_STATIC_LIB
         NAMES tbb
         HINTS
            ${CMAKE_CURRENT_BINARY_DIR}/tbb_release
            ${CMAKE_CURRENT_BINARY_DIR}/tbb_debug)
      if(TBB_STATIC_LIB)
         # create an imported target that also links the thread libraries via its interface
         set(TBB_BUILT_STATIC_LIB 1)
         find_package( Threads )
         list(APPEND linkinterface "${CMAKE_THREAD_LIBS_INIT}")
         list(APPEND linkinterface $<$<PLATFORM_ID:Linux>:rt>)     # Link "rt" library on Linux
         add_library(tbb STATIC IMPORTED)
         set_target_properties(tbb PROPERTIES
            IMPORTED_LOCATION "${TBB_STATIC_LIB}"
            INTERFACE_LINK_LIBRARIES "${linkinterface}")
         set(TBB_LIBRARY tbb)
      endif()
   endif()
endif()

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Red         "${Esc}[31m")
endif()

if(APPLE)
    message(STATUS "${Red}Please make sure that you have tbb 2020 installed, PaPILO < version 2 does not work with tbb 2021 and further. If the wrong version or none has been found, please specify TBB_DIR in you cmake call (i.e. cmake .. -DTBB_DIR=/path/to/tbb).${ColourReset}")
elseif(NOT TBB_LIBRARY)
    message(STATUS "${Red}If you have tbb 2020 installed and it has not been found, please specify TBB_DIR in you cmake call (i.e. cmake .. -DTBB_DIR=/path/to/tbb).${ColourReset}")
endif()

# TODO: modify FindTBB so that the version is saved in the corresponding variables
#if(TBB_LIBRARY)
#   file(READ "${TBB_LIBRARY}/tbb/tbb_stddef.h" _tbb_version_file)
#   string(REGEX REPLACE ".*#define TBB_VERSION_MAJOR ([0-9]+).*" "\\1"
#           TBB_VERSION_MAJOR "${_tbb_version_file}")
#   string(REGEX REPLACE ".*#define TBB_VERSION_MINOR ([0-9]+).*" "\\1"
#           TBB_VERSION_MINOR "${_tbb_version_file}")
#   string(REGEX REPLACE ".*#define TBB_INTERFACE_VERSION ([0-9]+).*" "\\1"
#           TBB_INTERFACE_VERSION "${_tbb_version_file}")
#   set(TBB_VERSION "${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR}")
#   message("${TBB_VERSION}")
#endif()
find_package_handle_standard_args(TBB REQUIRED_VARS TBB_LIBRARY)
