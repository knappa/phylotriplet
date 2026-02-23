add_executable(test_tree test_tree.cpp)

target_include_directories(test_tree PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}        # catch.hpp
)

target_link_libraries(test_tree PRIVATE phylotriplet)

find_package(Eigen3 REQUIRED)
target_link_libraries(test_tree PRIVATE Eigen3::Eigen)

add_test(NAME tree_tests COMMAND test_tree)