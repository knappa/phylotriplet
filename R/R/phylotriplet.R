#' Estimate phylogenetic tree from weighted trees file
#'
#' Reads a file containing Newick trees with weights and estimates a
#' phylogenetic tree. Returns a phylo object compatible with the APE package.
#'
#' @param filepath Path to file containing weighted trees
#'   (format: "newick_tree weight")
#'
#' @return An object of class "phylo" (APE format)
#'
#' @export
#' @importFrom ape read.tree
#'
#' @examples
#' \dontrun{
#' tree <- estimate_tree_from_file("weighted_trees.txt")
#' plot(tree)
#' write.tree(tree, "output.nwk")
#' }
estimate_tree_from_file <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }

  newick_string <- estimate_tree_from_file_rcpp(filepath)
  tree <- ape::read.tree(text = newick_string)

  attr(tree, "source_file") <- filepath
  attr(tree, "method") <- "triplet_spectral"

  return(tree)
}

#' Estimate phylogenetic tree from a string of weighted trees
#'
#' Parses a string containing Newick trees with weights and estimates a
#' phylogenetic tree. Returns a phylo object compatible with the APE package.
#'
#' @param content String containing weighted trees (format: "newick_tree weight",
#'   one per line)
#'
#' @return An object of class "phylo" (APE format)
#'
#' @export
#' @importFrom ape read.tree
#'
#' @examples
#' \dontrun{
#' tree <- estimate_tree_from_string("((A,B),C); 1.0\n((A,B),C); 2.0\n")
#' plot(tree)
#' }
estimate_tree_from_string <- function(content) {
  newick_string <- estimate_tree_from_string_rcpp(content)
  tree <- ape::read.tree(text = newick_string)

  attr(tree, "method") <- "triplet_spectral"

  return(tree)
}

#' Estimate phylogenetic tree from data.frame of Newick strings
#'
#' Estimates a phylogenetic tree from a data.frame containing Newick format
#' trees and their associated weights.
#'
#' @param trees_df Data frame with columns: newick (character), weight (numeric)
#'
#' @return An object of class "phylo" (APE format)
#'
#' @export
#' @importFrom ape read.tree
#'
#' @examples
#' \dontrun{
#' trees <- data.frame(
#'   newick = c("((human,chimp),gorilla);", "((human,chimp),orangutan);"),
#'   weight = c(0.95, 0.88),
#'   stringsAsFactors = FALSE
#' )
#' tree <- estimate_tree(trees)
#' plot(tree)
#' }
estimate_tree <- function(trees_df) {
  required_cols <- c("newick", "weight")
  if (!all(required_cols %in% names(trees_df))) {
    stop("trees_df must have columns: newick, weight")
  }

  if (nrow(trees_df) == 0) {
    stop("trees_df cannot be empty")
  }

  newick_string <- estimate_tree_rcpp(trees_df)
  tree <- ape::read.tree(text = newick_string)

  attr(tree, "input_tree_count") <- nrow(trees_df)
  attr(tree, "method") <- "triplet_spectral"

  return(tree)
}

#' Estimate tree from list of phylo objects with weights
#'
#' Estimates a phylogenetic tree from a list of phylo objects (or a data.frame
#' containing phylo objects) with associated weights.
#'
#' @param phylo_data Either a list where each element has 'tree' (phylo) and
#'   'weight', or a data.frame with 'tree' (list column of phylo objects)
#'   and 'weight' columns
#'
#' @return An object of class "phylo"
#'
#' @export
#' @importFrom ape write.tree read.tree
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree1 <- read.tree(text = "((human,chimp),gorilla);")
#' tree2 <- read.tree(text = "((human,chimp),orangutan);")
#'
#' phylo_list <- list(
#'   list(tree = tree1, weight = 0.95),
#'   list(tree = tree2, weight = 0.88)
#' )
#' result <- estimate_tree_from_phylo(phylo_list)
#' }
estimate_tree_from_phylo <- function(phylo_data) {
  if (is.data.frame(phylo_data)) {
    if (!all(c("tree", "weight") %in% names(phylo_data))) {
      stop("Data frame must have 'tree' and 'weight' columns")
    }
    trees <- phylo_data$tree
    weights <- phylo_data$weight
  } else if (is.list(phylo_data)) {
    trees <- lapply(phylo_data, function(x) {
      if ("tree" %in% names(x)) x$tree else x
    })
    weights <- sapply(phylo_data, function(x) {
      if ("weight" %in% names(x)) x$weight else 1.0
    })
  } else {
    stop("phylo_data must be a data.frame or list")
  }

  newick_strings <- sapply(trees, function(t) {
    if (!inherits(t, "phylo")) {
      stop("All trees must be phylo objects")
    }
    ape::write.tree(t)
  })

  trees_df <- data.frame(
    newick = newick_strings,
    weight = weights,
    stringsAsFactors = FALSE
  )

  estimate_tree(trees_df)
}

#' Read weighted trees file to data.frame
#'
#' Reads a file containing Newick trees with weights and returns a data.frame.
#'
#' @param filepath Path to weighted trees file (format: "newick_tree weight")
#'
#' @return Data frame with columns: newick (character), weight (numeric)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' trees <- read_weighted_trees("weighted_trees.txt")
#' head(trees)
#' summary(trees$weight)
#' }
read_weighted_trees <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  read_weighted_trees_rcpp(filepath)
}