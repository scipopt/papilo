
project(c_interface_testing LANGUAGES C)

add_executable(c_interface_testing "main.c")
#include_directories( "../../../../soplex/src")
include_directories( "../../../../papilo/src")
target_link_libraries(c_interface_testing fix-shared)
