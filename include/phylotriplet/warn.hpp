#pragma once
#include <iostream>

// Warning macro: defaults to std::cerr.
// Define PHYLOTRIPLET_RCPP_WARNINGS before including this header (and include
// Rcpp.h first) to route warnings through Rcpp::warning() instead.
// Define PHYLOTRIPLET_PY_WARNINGS before including this header (and include
// pybind11/pybind11.h first) to route warnings through Python's warnings
// module.
#if defined(PHYLOTRIPLET_RCPP_WARNINGS)
#define PHYLOTRIPLET_WARN(msg) Rcpp::warning(msg)
#elif defined(PHYLOTRIPLET_PY_WARNINGS)
#define PHYLOTRIPLET_WARN(msg) pybind11::module_::import("warnings").attr("warn")(msg)
#else
#define PHYLOTRIPLET_WARN(msg)  \
  do {                          \
    std::cerr << (msg) << '\n'; \
  } while (0)
#endif