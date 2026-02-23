# phylotriplet - R Package

Phylogenetic tree inference from weighted input trees.
All estimation functions return APE `phylo` objects for seamless integration
with the R phylogenetics ecosystem.

## Prerequisites

**System packages** (install before building):

```bash
# Debian/Ubuntu
apt install libeigen3-dev

# macOS (Homebrew)
brew install eigen

# Or point to a custom Eigen3 location
export EIGEN3_INCLUDE_DIR=/path/to/eigen3
```

**R packages:**

```r
install.packages(c("Rcpp", "ape"))
# Optional: devtools for development workflow
install.packages("devtools")
```

A C++14-capable compiler is also required (`gcc >= 5`, `clang >= 3.4`,
or MSVC 2017+).

## Installation

From the repository root:

```bash
cd R
Rscript -e "Rcpp::compileAttributes('.')"
R CMD INSTALL .
```

Or from within R (using devtools):

```r
setwd("R")
devtools::install()
```

## Development Workflow

```r
setwd("R")
devtools::load_all()      # load without installing
devtools::document()      # regenerate docs from roxygen2 comments
```

After any C++ change, regenerate the Rcpp wrappers before reloading:

```r
Rcpp::compileAttributes(".")
devtools::load_all()
```

## API

All functions return an APE `phylo` object unless noted.

### From a file

```r
tree <- estimate_tree_from_file("weighted_trees.txt")
plot(tree)
```

### From a string

```r
content <- "((human,chimp),gorilla); 0.95\n((human,chimp),orangutan); 0.88\n"
tree <- estimate_tree_from_string(content)
plot(tree)
```

### From a data frame

```r
trees <- data.frame(
  newick = c("((human,chimp),gorilla);", "((human,chimp),orangutan);"),
  weight = c(0.95, 0.88),
  stringsAsFactors = FALSE
)
tree <- estimate_tree(trees)
plot(tree)
```

### From a list of `phylo` objects

```r
library(ape)
tree1 <- read.tree(text = "((human,chimp),gorilla);")
tree2 <- read.tree(text = "((human,chimp),orangutan);")

result <- estimate_tree_from_phylo(list(
  list(tree = tree1, weight = 0.95),
  list(tree = tree2, weight = 0.88)
))
plot(result)
```

### Read weighted trees to a data frame

`read_weighted_trees()` parses a file without estimating a tree, returning a
data frame with columns `newick` (character) and `weight` (numeric).

```r
df <- read_weighted_trees("weighted_trees.txt")
head(df)
summary(df$weight)
```

## Input Format

One tree per line: `newick_tree weight`

```
((human,chimp),gorilla); 0.95
((human,chimp),orangutan); 0.88
((human,gorilla),orangutan); 0.82
```

Newick strings must end with `;`.

## Clean Rebuild

```bash
cd R
rm -f src/*.o src/*.so src/RcppExports.cpp R/RcppExports.R
Rscript -e "Rcpp::compileAttributes('.')"
R CMD INSTALL .
```

## See Also

- Main repository [README](../README.md) for algorithm details
- [APE package](https://cran.r-project.org/package=ape) for `phylo` object utilities