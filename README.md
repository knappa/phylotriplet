# phylotriplet

Phylogenetic tree estimation from weighted input trees using spectral graph decomposition.

This library provides a cross-language implementation (C++, Python, R) for estimating
phylogenetic trees from weighted input trees. It implements the BUILD (Aho et al.) algorithm for
inferring trees from rooted triplets, extended to handle weighted triplets.

## Input Format

The library accepts weighted trees in the format:

```
newick_tree weight
```

where each line contains a phylogenetic tree in Newick format, a space, and a numeric weight.

### Example: `weighted_trees.txt`

```
((A,B),C); 0.95
((A,B),D); 0.88
((A,C),D); 0.82
((E,F),G); 0.92
```

## Installation

See [QUICKSTART.md](QUICKSTART.md) for a quick guide and [BUILD.md](BUILD.md) for full
build and development instructions.

We use the [Eigen3](https://libeigen.gitlab.io/) library which must be installed separately:
```bash
apt install libeigen3-dev   # Debian/Ubuntu
brew install eigen          # macOS
```

### Python

```bash
cd python && pip install .
```

For development: `pip install -e ".[dev]"`

### R

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
trees = """
((A,B),C); 0.95
((A,B),D); 0.88
"""
newick = phylotriplet.estimate_tree_from_string(trees)

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

# Method 1: Direct from file (returns APE phylo object)
tree <- estimate_tree_from_file("weighted_trees.txt")
plot(tree)

# Method 2: From a string
trees = """
((A,B),C); 0.95
((A,B),D); 0.88
"""
tree <- estimate_tree_from_string(trees)

# Method 3: From a data frame of Newick strings
trees <- data.frame(
  newick = c("((A,B),C);", "((A,B),D);"),
  weight = c(0.95, 0.88),
  stringsAsFactors = FALSE
)
tree <- estimate_tree(trees)

# Method 4: From a list of phylo objects
tree1 <- read.tree(text = "((A,B),C);")
tree2 <- read.tree(text = "((A,B),D);")
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
auto trees = std::string("((A,B),C); 0.95\n((A,B),D); 0.88\n");
std::string tree = phylotriplet::estimate_tree_from_string(trees);

// From a vector of WeightedTree
std::vector<phylotriplet::WeightedTree> trees;
trees.push_back({"((A,B),C);", 0.95});
trees.push_back({"((A,B),D);", 0.88});
std::string newick = phylotriplet::estimate_tree(trees);
```

## APE Integration (R)

The R package is designed to work with APE

- Returns `phylo` objects directly
- Accepts both Newick strings and `phylo` objects as input
- Compatible with standard APE functions (`plot.phylo`, `write.tree`, `drop.tip`, etc.)
- Works with downstream packages (phytools, geiger, phangorn, etc.)

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

## Citation

If you use this software in your research, please cite it using the metadata in
[CITATION.cff](CITATION.cff). e.g.

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
