# phylotriplet

Phylogenetic tree estimation from weighted input trees using spectral triplet methods.

This library provides a cross-language implementation (C++, Python, R) for estimating
phylogenetic trees from weighted input trees. It implements the Aho et al. algorithm for
inferring trees from rooted triplets, extended to handle weighted triplets and to produce
polytomies when the input is not fully resolved. The spectral Laplacian (Fiedler vector)
is used to partition taxa recursively into a tree.

## Features

- **Core C++ implementation** for performance
- **Python bindings** via pybind11
- **R package** with seamless APE integration
- **File-based input** for large datasets
- **Flexible API** supporting both file and programmatic input

## Input Format

The library accepts weighted trees in a simple text format:

```
newick_tree weight
```

Each line contains a phylogenetic tree in Newick format, a space, and a numeric weight
(confidence, frequency, or posterior probability).

### Example: `weighted_trees.txt`

```
((human,chimp),gorilla); 0.95
((human,chimp),orangutan); 0.88
((human,gorilla),orangutan); 0.82
((mouse,rat),hamster); 0.92
```

Newick strings must end with `;`.

## Installation

See [QUICKSTART.md](QUICKSTART.md) for a quick guide and [BUILD.md](BUILD.md) for full
build and development instructions.

### Python

System prerequisite:
```bash
apt install libeigen3-dev   # Debian/Ubuntu
brew install eigen           # macOS
```

```bash
cd python && pip install .
```

For development: `pip install -e ".[dev]"`

### R

System prerequisite:
```bash
apt install libeigen3-dev   # Debian/Ubuntu
brew install eigen           # macOS
```

```bash
cd R
Rscript -e "Rcpp::compileAttributes('.')"
R CMD INSTALL .
```

## Usage

### Python

```python
import phylotriplet

# From a file
newick = phylotriplet.estimate_tree_from_file("weighted_trees.txt")
print(newick)

# From a string
content = "((human,chimp),gorilla); 0.95\n((human,chimp),orangutan); 0.88\n"
newick = phylotriplet.estimate_tree_from_string(content)

# Read trees without estimating
trees = phylotriplet.read_weighted_trees_from_file("weighted_trees.txt")
print(f"Loaded {len(trees)} trees")

# Visualise with ete3
from ete3 import Tree
tree = Tree(newick)
print(tree)
```

### R

```r
library(phylotriplet)
library(ape)

# Method 1: Direct from file (returns phylo object)
tree <- estimate_tree_from_file("weighted_trees.txt")
plot(tree)

# Method 2: From a string
content <- "((human,chimp),gorilla); 0.95\n((human,chimp),orangutan); 0.88\n"
tree <- estimate_tree_from_string(content)

# Method 3: From a data frame of Newick strings
trees <- data.frame(
  newick = c("((human,chimp),gorilla);", "((human,chimp),orangutan);"),
  weight = c(0.95, 0.88),
  stringsAsFactors = FALSE
)
tree <- estimate_tree(trees)

# Method 4: From a list of phylo objects
tree1 <- read.tree(text = "((human,chimp),gorilla);")
tree2 <- read.tree(text = "((human,chimp),orangutan);")
tree <- estimate_tree_from_phylo(list(
  list(tree = tree1, weight = 0.95),
  list(tree = tree2, weight = 0.88)
))

# Result works with all APE functions
write.tree(tree, "result.nwk")
is.rooted(tree)
Ntip(tree)
```

### C++

```cpp
#include "phylotriplet/triplet_estimator.hpp"

// From a file
std::string newick = phylotriplet::estimate_tree_from_file("weighted_trees.txt");

// From a string
std::string newick = phylotriplet::estimate_tree_from_string(content);

// From a vector of WeightedTree
std::vector<phylotriplet::WeightedTree> trees;
trees.push_back({"((human,chimp),gorilla);", 0.95});
trees.push_back({"((human,chimp),orangutan);", 0.88});
std::string newick = phylotriplet::estimate_tree(trees);
```

## APE Integration (R)

The R package is designed to work seamlessly with APE (Analysis of Phylogenetics and Evolution):

- Returns `phylo` objects directly
- Accepts both Newick strings and `phylo` objects as input
- Compatible with all APE functions (`plot.phylo`, `write.tree`, `drop.tip`, etc.)
- Works with downstream packages (phytools, geiger, phangorn, etc.)

## Project Structure

```
phylotriplet/
├── include/phylotriplet/
│   ├── triplet_estimator.hpp  # Top-level API: estimate_tree*
│   ├── gen_tree.hpp           # WeightedTriplet, spectral recursive tree building
│   ├── spectral_partition.hpp # Fiedler-vector Laplacian partitioning (Eigen3)
│   ├── tree.hpp               # TreeNode / PhyloTree (Newick parse/emit, triplets)
│   └── file_io.hpp            # File parser, WeightedTree struct
├── tools/
│   └── newick2triplets.cpp    # CLI: enumerate triplets from a Newick tree
├── python/
│   ├── bindings.cpp           # pybind11 bindings
│   └── setup.py               # Python package setup
├── R/
│   ├── DESCRIPTION            # R package metadata
│   ├── NAMESPACE              # R exports
│   ├── R/phylotriplet.R       # R wrapper functions (returns phylo objects)
│   └── src/rcpp_bindings.cpp  # Rcpp bindings
├── tests/
│   └── test_tree.cpp          # C++ unit tests (Catch2)
└── examples/
    └── weighted_trees.txt     # Example data
```

## Development

### Python

```bash
cd python
pip install -e . --force-reinstall --no-deps   # rebuild after C++ changes
```

### R

```r
Rcpp::compileAttributes(".")   # after C++ changes
devtools::load_all()
devtools::test()
```

### C++ tests

```bash
cmake -B build -DBUILD_TESTS=ON && cmake --build build && ctest --test-dir build
```

### Command-line tool

```bash
cmake -B build -DBUILD_TOOLS=ON && cmake --build build --target newick2triplets
build/newick2triplets "((a,b),c);"              # from arg
build/newick2triplets --file examples/weighted_trees.txt  # from file
```

## Requirements

### C++
- C++14 compatible compiler
- Eigen3 (`apt install libeigen3-dev` / `brew install eigen`)
- CMake 3.14+ (optional, for tests/tools only)

### Python
- Python 3.7+
- pybind11 (fetched automatically as a build dependency)
- Eigen3 (system package, see above)

### R
- R 3.5+
- Rcpp (>= 1.0.0)
- ape (>= 5.0)
- Eigen3 (system package, see above)

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use this software in your research, please cite it using the metadata in
[CITATION.cff](CITATION.cff). A BibTeX entry:

```bibtex
@software{knapp_phylotriplet,
  author    = {Knapp, Adam C. and Chifman, Julia},
  title     = {phylotriplet},
  url       = {https://github.com/knappa/phylotriplet},
  license   = {MIT},
}
```

## Authors

- Adam C. Knapp (adam.knapp@medicine.ufl.edu)
- Julia Chifman (chifman@american.edu)

## Contributing

Contributions are welcome. Please fork the repository, create a feature branch,
add tests, and submit a pull request.

Issues: https://github.com/knappa/phylotriplet/issues