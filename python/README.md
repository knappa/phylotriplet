# phylotriplet - Python Package

Phylogenetic tree inference from weighted input trees.

## Prerequisites

Eigen3 must be installed at the system level before building:

```bash
# Debian/Ubuntu
apt install libeigen3-dev

# macOS (Homebrew)
brew install eigen

# Or set EIGEN3_INCLUDE_DIR to a custom path
export EIGEN3_INCLUDE_DIR=/path/to/eigen3
```

## Installation

```bash
pip install .
```

`pybind11` is fetched automatically as a build dependency — no separate install needed.

For development (includes `pytest`, `ete3`):

```bash
pip install -e ".[dev]"
```

## API

### Tree estimation

```python
import phylotriplet

# From a weighted-trees file
newick = phylotriplet.estimate_tree_from_file("weighted_trees.txt")

# From a weighted-trees string
content = "((human,chimp),gorilla); 0.95\n((human,chimp),orangutan); 0.88\n"
newick = phylotriplet.estimate_tree_from_string(content)
```

Both return a Newick string.

### I/O utilities

```python
# Read a weighted-trees file into a list of WeightedTree objects
trees = phylotriplet.read_weighted_trees_from_file("weighted_trees.txt")

# Parse a weighted-trees string into a list of WeightedTree objects
trees = phylotriplet.read_weighted_trees_from_string(content)
```

## Input Format

One tree per line: `newick_tree weight`

```
((human,chimp),gorilla); 0.95
((human,chimp),orangutan); 0.88
((human,gorilla),orangutan); 0.82
```

Newick strings must end with `;`.

## Quick Start

```python
import phylotriplet

newick = phylotriplet.estimate_tree_from_file("examples/weighted_trees.txt")
print(newick)

# Visualise with ete3
from ete3 import Tree
tree = Tree(newick)
print(tree)
```

See the main [README.md](../README.md) for more usage examples.