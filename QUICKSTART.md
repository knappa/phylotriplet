# Quick Start Guide

## Prerequisites

Install Eigen3 before building either package:

```bash
apt install libeigen3-dev   # Debian/Ubuntu
brew install eigen           # macOS
# custom path: export EIGEN3_INCLUDE_DIR=/path/to/eigen3
```

## Installation

### Python
```bash
cd python
pip install .
```

### R
```bash
cd R
Rscript -e "Rcpp::compileAttributes('.')"
R CMD INSTALL .
```

## Try It

```bash
python examples/test_python.py
Rscript examples/test_r.R
```

## Use Your Own Data

Create `my_trees.txt` (one Newick tree + weight per line):

```
((human,chimp),gorilla); 0.95
((human,chimp),orangutan); 0.88
((human,gorilla),orangutan); 0.82
((mouse,rat),hamster); 0.92
```

**Python:**
```python
import phylotriplet as pt

newick = pt.estimate_tree_from_file("my_trees.txt")
print(newick)

# Visualise with ete3
from ete3 import Tree
print(Tree(newick))
```

**R:**
```r
library(phylotriplet)

tree <- estimate_tree_from_file("my_trees.txt")
plot(tree)
write.tree(tree, "result.nwk")
```

## Input Format

Each line: `newick_tree weight`

- Newick string must end with `;`
- Weight is the last whitespace-delimited token on the line
- Blank lines and lines starting with `#` are ignored

```
((A,B),C); 0.95
((A:0.1,B:0.2):0.3,C:0.4); 0.75
```

## API Summary

### Python

| Function | Input | Returns |
|---|---|---|
| `estimate_tree_from_file(filepath)` | path to weighted-trees file | Newick string |
| `estimate_tree_from_string(content)` | weighted-trees string | Newick string |
| `read_weighted_trees_from_file(filepath)` | path to weighted-trees file | list of `WeightedTree` |
| `read_weighted_trees_from_string(content)` | weighted-trees string | list of `WeightedTree` |

### R

| Function | Input | Returns |
|---|---|---|
| `estimate_tree_from_file(filepath)` | path to weighted-trees file | `phylo` |
| `estimate_tree_from_string(content)` | weighted-trees string | `phylo` |
| `estimate_tree(trees_df)` | data.frame with `newick`, `weight` cols | `phylo` |
| `estimate_tree_from_phylo(phylo_data)` | list of `list(tree=phylo, weight=...)` | `phylo` |
| `read_weighted_trees(filepath)` | path to weighted-trees file | data.frame |

## Troubleshooting

**Eigen3 not found:**
```bash
apt install libeigen3-dev   # Linux
brew install eigen           # macOS
export EIGEN3_INCLUDE_DIR=/path/to/eigen3  # custom
```

**Python: module not found** — run `pip install .` from the `python/` directory.

**R: package not installed** — run `R CMD INSTALL .` from the `R/` directory.

**Rebuild after C++ changes:**
```bash
# Python
cd python && pip install -e . --force-reinstall --no-deps

# R
Rscript -e "Rcpp::compileAttributes('.')" && R CMD INSTALL R/
```

For full build and development instructions see [BUILD.md](BUILD.md).

## Issues

https://github.com/knappa/phylotriplet/issues