# Build and Development Guide

Detailed instructions for building, testing, and developing phylotriplet.

## Prerequisites

### All platforms
- C++14 compatible compiler (GCC 5+, Clang 3.4+, MSVC 2015+)
- **Eigen3** — required by the spectral partitioning algorithm

```bash
apt install libeigen3-dev   # Debian/Ubuntu
brew install eigen           # macOS (Homebrew)
# Or set EIGEN3_INCLUDE_DIR=/path/to/eigen3
```

### Python
- Python 3.7+
- pip (pybind11 is fetched automatically as a build dependency)

### R
- R 3.5+
- Rcpp, ape (`install.packages(c("Rcpp", "ape"))`)
- devtools (optional, for development workflow)

## Building

### Python Package

```bash
cd python
pip install .
```

Development (editable) install:
```bash
cd python
pip install -e ".[dev]"
```

After modifying C++ headers:
```bash
pip install -e . --force-reinstall --no-deps
```

Verify:
```bash
python -c "import phylotriplet; print(phylotriplet.__doc__)"
```

### R Package

```bash
cd R
Rscript -e "Rcpp::compileAttributes('.')"
R CMD INSTALL .
```

From the R console (iterative development):
```r
setwd("R")
devtools::install()      # full install
# or
devtools::load_all()     # load without installing
```

Verify:
```r
library(phylotriplet)
?estimate_tree_from_file
```

### C++ Tests

```bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
ctest --test-dir build
```

### Command-line Tool (`newick2triplets`)

```bash
cmake -B build -DBUILD_TOOLS=ON
cmake --build build --target newick2triplets
```

Usage:
```bash
build/newick2triplets "((a,b),c);"               # single Newick from arg (weight=1)
echo "((a,b),c);" | build/newick2triplets         # from stdin
build/newick2triplets --file examples/weighted_trees.txt  # weighted-trees file
```

## Development Workflow

### Python

1. Edit C++ headers in `include/phylotriplet/`
2. Update `python/bindings.cpp` if the API changes
3. Rebuild: `pip install -e . --force-reinstall --no-deps`
4. Test: `python examples/test_python.py` or `pytest tests/`

### R

1. Edit C++ headers in `include/phylotriplet/`
2. Update `R/src/rcpp_bindings.cpp` if the API changes
3. Regenerate Rcpp exports: `Rcpp::compileAttributes(".")`
4. Update docs: `devtools::document()`
5. Reload: `devtools::load_all()`
6. Test: `source("examples/test_r.R")` or `devtools::test()`

## Testing

### Python Tests
```bash
pytest tests/test_python.py
# or run the example script
python examples/test_python.py
```

### R Tests
```r
devtools::test()
# or run the example script
source("examples/test_r.R")
```

### C++ Tests
```bash
ctest --test-dir build --output-on-failure
```

## Package Checks

### Python
```bash
cd python
pip install build twine
python -m build
twine check dist/*
```

### R
```bash
cd R
R CMD check .
# or
devtools::check()
```

## Documentation

### Python
Docstrings are written in `python/bindings.cpp` as the second argument to each `m.def()` call.

### R
Roxygen2 comments in `R/R/phylotriplet.R`. After editing:
```r
devtools::document()
```

## Clean Rebuild

```bash
# Python
rm -rf python/build python/*.egg-info
cd python && pip install -e . --force-reinstall

# R
rm -f R/src/*.o R/src/*.so
cd R && R CMD INSTALL .
```

## Distribution

### Python (PyPI)
```bash
cd python
python -m build
twine upload dist/*
```

### R (CRAN)
1. `R CMD check --as-cran .`
2. Fix all warnings and notes
3. Update version in `R/DESCRIPTION`
4. Submit to CRAN

## Common Issues

### "Eigen3 not found"
```bash
apt install libeigen3-dev        # Linux
brew install eigen               # macOS
export EIGEN3_INCLUDE_DIR=...    # custom path
```

### Python: "fatal error: pybind11/pybind11.h: No such file"
pybind11 should be fetched automatically. If not, check that pip >= 21.0:
```bash
pip install --upgrade pip && pip install .
```

### R: "function 'Rcpp_precious_remove' not provided by package 'Rcpp'"
```r
update.packages("Rcpp")
Rcpp::compileAttributes(".")
```

### Linking errors after modifying C++ code
```bash
# Python
rm -rf python/build python/*.egg-info
cd python && pip install -e . --force-reinstall

# R
rm -f R/src/*.o R/src/*.so
cd R && R CMD INSTALL .
```

## Support

Issues: https://github.com/knappa/phylotriplet/issues