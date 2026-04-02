cmake_minimum_required(VERSION 3.14)
project(phylotriplet VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(BUILD_PYTHON "Build Python bindings" OFF)
option(BUILD_TESTS "Build C++ tests" OFF)
option(BUILD_TOOLS "Build command-line tools" OFF)

# Header-only library
add_library(phylotriplet INTERFACE)
target_include_directories(phylotriplet INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/LBFGSpp>
    $<INSTALL_INTERFACE:include>
)

# Python bindings (handled by setup.py, this is mainly for IDE support)
if(BUILD_PYTHON)
    message(STATUS "Python bindings should be built using setup.py")
    message(STATUS "  cd python && pip install .")
endif()

# C++ tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Command-line tools
if(BUILD_TOOLS)
    find_package(Eigen3 REQUIRED)
    add_executable(newick2triplets tools/newick2triplets.cpp)
    target_link_libraries(newick2triplets PRIVATE phylotriplet Eigen3::Eigen)
    install(TARGETS newick2triplets RUNTIME DESTINATION bin)
endif()

# Installation
install(DIRECTORY include/phylotriplet
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")

# Package config
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/phylotripletConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/phylotripletConfigVersion.cmake"
    DESTINATION lib/cmake/phylotriplet
)
