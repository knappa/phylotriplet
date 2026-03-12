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

### Sequence alignment input (Python)

```python
import phylotriplet as pt

# Read a NEXUS file (single or multi-gene)
aln = pt.read_nexus_file("mydata.nex")
print(aln)                   # MultiGeneAlignment(num_genes=2, num_taxa=10)
print(aln.all_taxon_names()) # sorted list of taxon names

# Read one or more PHYLIP files
aln = pt.read_phylip_files(["gene1.phy", "gene2.phy"],
                            gene_names=["COI", "CYTB"])

# Compute aggregate triplet site-pattern counts
patterns = aln.compute_site_patterns()   # list of TuplePatterns
for p in patterns:
    # p.taxa       — sorted list of 3 taxon names
    # p.n          — list of 64 floats; p.n[i*16 + j*4 + k] counts sites
    #                where taxa[0]=i, taxa[1]=j, taxa[2]=k (0=A,1=C,2=G,3=T)
    # p.n_total    — total alignment columns
    # p.n_counted  — columns with no gap in any of the three taxa
    print(p.taxa, p.n_counted)
```

### Sequence alignment input (R)

```r
library(phylotriplet)

# Read a NEXUS file
aln <- read_nexus("mydata.nex")
print(aln)
# PhylotripletAlignment
#   2 gene(s), 10 taxon/taxa
#   Taxa: Aardvark, Bear, Cat, ...

# Read one or more PHYLIP files
aln <- read_phylip_files(c("gene1.phy", "gene2.phy"),
                          gene_names = c("COI", "CYTB"))

# Compute aggregate triplet site-pattern counts
pats <- compute_site_patterns(aln)
# Returns a data frame with columns:
#   taxon1, taxon2, taxon3  — taxon names (lexicographic order)
#   n_total, n_counted      — alignment column counts
#   n                       — list column; each element is numeric(64)

head(pats[, c("taxon1", "taxon2", "taxon3", "n_total", "n_counted")])

# Access the 64-vector for row 1:
#   pats$n[[1]][a*16 + b*4 + c + 1]  (1-indexed; 0=A,1=C,2=G,3=T)
pats$n[[1]]

# Fraction of non-gap sites per triplet:
pats$n_counted / pats$n_total
```

## APE Integration (R)

The R package is designed to work with APE

- Returns `phylo` objects directly
- Accepts both Newick strings and `phylo` objects as input
- Compatible with standard APE functions (`plot.phylo`, `write.tree`, `drop.tip`, etc.)
- Works with downstream packages (phytools, geiger, phangorn, etc.)

## Testing

### C++

```bash
cmake -B build -DBUILD_TESTS=ON && cmake --build build && ctest --test-dir build
```

Two test suites run under [Catch2](https://github.com/catchorg/Catch2):

| Executable | Source | What it covers |
|---|---|---|
| `test_tree` | `tests/test_tree.cpp` | Newick parse/emit, triplets, Robinson-Foulds, prune, LCA, traversal |
| `test_sequence` | `tests/test_sequence.cpp` | `char_to_nucleotide`, `nucleotide_distribution`, PHYLIP/NEXUS parsing, site-pattern computation |

To run a single suite:

```bash
./build/tests/test_tree
./build/tests/test_sequence
```

### Python

Requires pytest (`pip install pytest`).

```bash
# From the repo root:
pytest tests/test_bindings_python.py -v
```

`tests/test_bindings_python.py` covers the Python bindings for sequence alignment:
`read_nexus_file`, `read_phylip_file`, `read_phylip_files`,
`MultiGeneAlignment`, and `TuplePatterns`.

### R

```r
# From inside the R/ directory (or any R session with the package loaded):
devtools::test("R/")
```

Testthat tests live in `R/tests/testthat/test-sequence.R` and cover
`read_nexus`, `read_phylip_files`, `print.PhylotripletAlignment`, and
`compute_site_patterns`.

Alternatively, install and run from the shell:

```bash
cd R && R CMD check . --no-manual
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
