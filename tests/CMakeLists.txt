find_package(Eigen3 REQUIRED)

add_executable(test_tree test_tree.cpp)
target_include_directories(test_tree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_tree PRIVATE phylotriplet Eigen3::Eigen)
add_test(NAME tree_tests COMMAND test_tree)

add_executable(test_sequence test_sequence.cpp)
target_include_directories(test_sequence PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_sequence PRIVATE phylotriplet Eigen3::Eigen)
add_test(NAME sequence_tests COMMAND test_sequence)