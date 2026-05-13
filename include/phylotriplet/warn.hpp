#pragma once
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <string>

// Printf-style string formatter for use with PHYLOTRIPLET_WARN.
// Avoids std::format (C++20) and keeps C++14 compatibility.
inline std::string phylotriplet_warn_format(const char* fmt, ...) {
  char buf[256];
  va_list args;
  va_start(args, fmt);
  vsnprintf(buf, sizeof(buf), fmt, args);
  va_end(args);
  return std::string(buf);
}

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