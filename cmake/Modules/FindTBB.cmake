include(FindPackageHandleStandardArgs)

# Define search paths based on user input and environment variables
set(TBB_SEARCH_DIR ${TBB_LIBRARY_DIR} ${TBB_ROOT_DIR} $ENV{TBB_INSTALL_DIR} $ENV{TBBROOT})

# for windows add additional default search paths
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
   set(TBB_DEFAULT_SEARCH_DIR  "C:/Program Files/Intel/TBB"
                               "C:/Program Files (x86)/Intel/TBB")
endif()

# try to find the library
find_library(TBB_LIBRARY
            NAMES tbb libtbb.so.2
            HINTS ${TBB_SEARCH_DIR}
            PATHS ${TBB_DEFAULT_SEARCH_DIR})

find_package_handle_standard_args(TBB REQUIRED_VARS TBB_LIBRARY)
