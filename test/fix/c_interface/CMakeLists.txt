project(c_interface_testing LANGUAGES C)

add_executable(c_interface_testing "main.c")
#include_directories( "../../../../soplex/src")
include_directories("${PROJECT_SOURCE_DIR}/src")
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS program_options iostreams REQUIRED)
target_link_libraries(c_interface_testing fix-shared papilo ${PAPILOLIB_TARGET} Boost::iostreams)
