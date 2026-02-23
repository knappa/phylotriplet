# Installation Fix Guide

If you got a `ModuleNotFoundError: No module named 'pybind11'` error, follow these updated instructions.

## Python Installation (Fixed)

The project now includes a `pyproject.toml` file that automatically handles pybind11 installation.

### Option 1: Simple Installation (Recommended)

```bash
cd python
pip install .
```

That's it! The `pyproject.toml` file tells pip to install pybind11 first, then build the package.

### Option 2: Development Mode

```bash
cd python
pip install -e .
```

### Option 3: If You Still Have Issues

If you're using an older version of pip, upgrade it first:

```bash
pip install --upgrade pip
cd python
pip install .
```

### Verify Installation

```bash
python -c "import phylotriplet; print('Success!')"
```

### Test It

```bash
cd ../examples
python test_python.py
```

## What Changed?

The project now includes `python/pyproject.toml` which specifies build dependencies:

```toml
[build-system]
requires = ["setuptools>=42", "wheel", "pybind11>=2.6.0"]
build-backend = "setuptools.build_meta"
```

This tells pip: "Before trying to build this package, install pybind11 first."

## R Installation (No Changes)

R installation works as before:

```bash
cd R
Rscript -e "Rcpp::compileAttributes('.')"
R CMD INSTALL .
```

Or from R:
```r
devtools::document()
devtools::install()
```

## Need More Help?

If you still encounter issues:

1. **Check Python version**: `python --version` (need 3.7+)
2. **Check pip version**: `pip --version` (upgrade if older than 21.0)
3. **Try in a fresh virtual environment**:
   ```bash
   python -m venv myenv
   source myenv/bin/activate  # On Windows: myenv\Scripts\activate
   cd python
   pip install .
   ```

4. **Check compiler**: Make sure you have a C++14 compatible compiler
   - Linux: `gcc --version` (need 5.0+)
   - macOS: Xcode command line tools
   - Windows: Visual Studio 2015+

## Quick Test After Installation

### Python
```python
import phylotriplet

# Estimate from a string
content = "((A,B),C); 0.95\n((A,C),B); 0.40\n"
result = phylotriplet.estimate_tree_from_string(content)
print(f"Success! Tree: {result}")
```

### R
```r
library(phylotriplet)

# Create a simple test
trees <- data.frame(
  newick = "((A,B),C);",
  weight = 0.95,
  stringsAsFactors = FALSE
)
result <- estimate_tree(trees)
print(paste("Success! Tree has", ape::Ntip(result), "tips"))
```

Both should run without errors!
