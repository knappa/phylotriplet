find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(phylotriplet bindings.cpp)

target_link_libraries(phylotriplet PRIVATE phylotriplet)

# Installation
install(TARGETS phylotriplet
  LIBRARY DESTINATION .
  RUNTIME DESTINATION .)
