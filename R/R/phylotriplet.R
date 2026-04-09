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

# -- Sequence alignment ---------------------------------------------------------

#' Read a NEXUS alignment file
#'
#' Parses a NEXUS file and returns a PhylotripletAlignment object. Supports
#' multiple DATA/CHARACTERS blocks (one per gene), a single DATA block split
#' by a SETS block with CHARSET commands, sequential and interleaved MATRIX
#' formats, and DNA/RNA datatypes.
#'
#' @param filepath Path to the NEXUS file
#'
#' @return A PhylotripletAlignment object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln <- read_nexus("mydata.nex")
#' print(aln)
#' pats <- triplet_site_patterns(aln)
#' head(pats[, c("taxon1", "taxon2", "taxon3", "n_total", "n_counted")])
#' }
read_nexus <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  ptr <- read_nexus_rcpp(filepath)
  structure(list(.ptr = ptr), class = "PhylotripletAlignment")
}

#' Read one or more PHYLIP alignment files
#'
#' Reads one or more PHYLIP-format files as a multi-gene alignment and returns
#' a PhylotripletAlignment object.
#'
#' @param filepaths Character vector of file paths, one per gene
#' @param gene_names Optional character vector of gene names; defaults to
#'   filename stems
#' @param interleaved Logical; set TRUE for interleaved PHYLIP format
#'   (applied uniformly to all files)
#'
#' @return A PhylotripletAlignment object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln <- read_phylip_files(c("gene1.phy", "gene2.phy"))
#' pats <- triplet_site_patterns(aln)
#' }
read_phylip_files <- function(filepaths, gene_names = character(0),
                               interleaved = FALSE) {
  missing <- filepaths[!file.exists(filepaths)]
  if (length(missing) > 0) {
    stop("File(s) not found: ", paste(missing, collapse = ", "))
  }
  ptr <- read_phylip_files_rcpp(filepaths, gene_names, interleaved)
  structure(list(.ptr = ptr), class = "PhylotripletAlignment")
}

#' Print a PhylotripletAlignment
#'
#' @param x A PhylotripletAlignment object
#' @param ... Ignored
#'
#' @export
print.PhylotripletAlignment <- function(x, ...) {
  ng <- alignment_num_genes_rcpp(x$.ptr)
  taxa <- alignment_taxon_names_rcpp(x$.ptr)
  cat("PhylotripletAlignment\n")
  cat(" ", ng, "gene(s),", length(taxa), "taxon/taxa\n")
  if (length(taxa) > 0) {
    preview <- if (length(taxa) <= 8) taxa else c(taxa[1:8], "...")
    cat("  Taxa:", paste(preview, collapse = ", "), "\n")
  }
  invisible(x)
}

# -- Gene access ---------------------------------------------------------------

#' Get gene names from an alignment
#'
#' @param alignment A PhylotripletAlignment object.
#' @param ... Ignored.
#'
#' @return Character vector of gene names in order.
#'
#' @export
gene_names <- function(alignment, ...) UseMethod("gene_names")

#' @export
gene_names.PhylotripletAlignment <- function(alignment, ...) {
  alignment_gene_names_rcpp(alignment$.ptr)
}

#' Extract a single gene from an alignment
#'
#' Returns a PhylotripletGene object that can be passed to
#' \code{\link{triplet_site_patterns}}, \code{\link{triplet_site_pattern}},
#' \code{\link{quartet_site_patterns}}, \code{\link{quartet_site_pattern}},
#' and \code{\link{site_pattern}}.
#'
#' @param alignment A PhylotripletAlignment object.
#' @param idx Either a 1-based integer index or a character gene name.
#' @param ... Ignored.
#'
#' @return A PhylotripletGene object.
#'
#' @export
get_gene <- function(alignment, idx, ...) UseMethod("get_gene")

#' @export
get_gene.PhylotripletAlignment <- function(alignment, idx, ...) {
  ptr <- get_gene_rcpp(alignment$.ptr, idx)
  structure(list(.ptr = ptr), class = "PhylotripletGene")
}

#' Print a PhylotripletGene
#'
#' @param x A PhylotripletGene object
#' @param ... Ignored
#'
#' @export
print.PhylotripletGene <- function(x, ...) {
  nm <- gene_name_rcpp(x$.ptr)
  taxa <- gene_taxon_names_rcpp(x$.ptr)
  cat("PhylotripletGene\n")
  cat("  Gene:", nm, "\n")
  cat(" ", length(taxa), "taxon/taxa\n")
  if (length(taxa) > 0) {
    preview <- if (length(taxa) <= 8) taxa else c(taxa[1:8], "...")
    cat("  Taxa:", paste(preview, collapse = ", "), "\n")
  }
  invisible(x)
}

#' Compute triplet site-pattern counts from an alignment
#'
#' Computes aggregate site-pattern counts across all genes for every n choose 3
#' taxon triplet.  Contributions from a gene are included only for genes that
#' contain all three taxa; n_total and n_counted are summed accordingly.
#'
#' @param alignment A PhylotripletAlignment object returned by
#'   \code{\link{read_nexus}} or \code{\link{read_phylip_files}}
#' @param ... Ignored
#'
#' @return A data frame with one row per taxon triplet and columns:
#'   \describe{
#'     \item{taxon1, taxon2, taxon3}{Taxon names (lexicographic order)}
#'     \item{n_total}{Total alignment columns across contributing genes}
#'     \item{n_counted}{Columns with no gap in any of the three taxa}
#'     \item{site_pattern}{List column; each element is a numeric vector of length 64.
#'       \code{site_pattern[[i]][a*16 + b*4 + c + 1]} (1-indexed) is the
#'       fractional count of sites where taxon1 shows base a, taxon2 shows b,
#'       taxon3 shows c, with 0=A, 1=C, 2=G, 3=T.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln  <- read_nexus("mydata.nex")
#' pats <- triplet_site_patterns(aln)
#'
#' # Access the 64-vector for the first triplet
#' pats$site_pattern[[1]]
#'
#' # Fraction of sites with no gap (per triplet)
#' pats$n_counted / pats$n_total
#' }
triplet_site_patterns <- function(alignment, ...) {
  UseMethod("triplet_site_patterns")
}

#' @export
triplet_site_patterns.PhylotripletAlignment <- function(alignment, ...) {
  raw <- triplet_site_patterns_rcpp(alignment$.ptr)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- lapply(seq_len(nrow(df)), function(i) raw$site_pattern[i, ])
  df
}

#' @export
triplet_site_patterns.PhylotripletGene <- function(alignment, ...) {
  raw <- gene_triplet_site_patterns_rcpp(alignment$.ptr)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- lapply(seq_len(nrow(df)), function(i) raw$site_pattern[i, ])
  df
}

#' Compute site-pattern counts for a specific taxon triplet
#'
#' Like \code{\link{triplet_site_patterns}} but restricted to a single
#' user-specified triplet of taxa.  Genes that do not contain all three taxa
#' are skipped with a warning (multi-gene alignment) or an error is raised
#' (single gene).  The taxon names may be given in any order.
#'
#' @param alignment A PhylotripletAlignment or PhylotripletGene object.
#' @param taxon1,taxon2,taxon3 Character scalars naming the three taxa.
#' @param ... Ignored.
#'
#' @return A one-row data frame with the same columns as
#'   \code{\link{triplet_site_patterns}}: \code{taxon1}, \code{taxon2},
#'   \code{taxon3}, \code{n_total}, \code{n_counted}, and \code{site_pattern}
#'   (a list column containing one numeric vector of length 64).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln  <- read_nexus("mydata.nex")
#' pats <- triplet_site_pattern(aln, "human", "chimp", "gorilla")
#' pats$site_pattern[[1]]
#' }
triplet_site_pattern <- function(alignment, taxon1, taxon2, taxon3, ...) {
  UseMethod("triplet_site_pattern")
}

#' @export
triplet_site_pattern.PhylotripletAlignment <- function(alignment, taxon1,
                                                       taxon2, taxon3, ...) {
  raw <- triplet_site_pattern_rcpp(alignment$.ptr, taxon1, taxon2, taxon3)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- list(raw$site_pattern)
  df
}

#' @export
triplet_site_pattern.PhylotripletGene <- function(alignment, taxon1,
                                                  taxon2, taxon3, ...) {
  raw <- gene_triplet_site_pattern_rcpp(alignment$.ptr, taxon1, taxon2, taxon3)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- list(raw$site_pattern)
  df
}

#' Compute quartet site-pattern counts from an alignment
#'
#' Computes aggregate site-pattern counts across all genes for every n choose 4
#' taxon quartet.  Contributions from a gene are included only for genes that
#' contain all four taxa; n_total and n_counted are summed accordingly.
#'
#' @param alignment A PhylotripletAlignment object returned by
#'   \code{\link{read_nexus}} or \code{\link{read_phylip_files}}
#' @param ... Ignored
#'
#' @return A data frame with one row per taxon quartet and columns:
#'   \describe{
#'     \item{taxon1, taxon2, taxon3, taxon4}{Taxon names (lexicographic order)}
#'     \item{n_total}{Total alignment columns across contributing genes}
#'     \item{n_counted}{Columns with no gap in any of the four taxa}
#'     \item{site_pattern}{List column; each element is a numeric vector of length 256.
#'       \code{site_pattern[[i]][a*64 + b*16 + c*4 + d + 1]} (1-indexed) is the
#'       fractional count of sites where taxon1=a, taxon2=b, taxon3=c, taxon4=d,
#'       with 0=A, 1=C, 2=G, 3=T.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln  <- read_nexus("mydata.nex")
#' pats <- quartet_site_patterns(aln)
#' pats$site_pattern[[1]]
#' }
quartet_site_patterns <- function(alignment, ...) {
  UseMethod("quartet_site_patterns")
}

#' @export
quartet_site_patterns.PhylotripletAlignment <- function(alignment, ...) {
  raw <- quartet_site_patterns_rcpp(alignment$.ptr)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    taxon4    = raw$taxon4,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- lapply(seq_len(nrow(df)), function(i) raw$site_pattern[i, ])
  df
}

#' @export
quartet_site_patterns.PhylotripletGene <- function(alignment, ...) {
  raw <- gene_quartet_site_patterns_rcpp(alignment$.ptr)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    taxon4    = raw$taxon4,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- lapply(seq_len(nrow(df)), function(i) raw$site_pattern[i, ])
  df
}

#' Compute site-pattern counts for a specific taxon quartet
#'
#' Like \code{\link{quartet_site_patterns}} but restricted to a single
#' user-specified quartet of taxa.  Genes that do not contain all four taxa
#' are skipped with a warning.  The taxon names may be given in any order.
#'
#' @param alignment A PhylotripletAlignment object.
#' @param taxon1,taxon2,taxon3,taxon4 Character scalars naming the four taxa.
#' @param ... Ignored.
#'
#' @return A one-row data frame with the same columns as
#'   \code{\link{quartet_site_patterns}}: \code{taxon1}, \code{taxon2},
#'   \code{taxon3}, \code{taxon4}, \code{n_total}, \code{n_counted}, and
#'   \code{site_pattern} (a list column containing one numeric vector of length 256).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln  <- read_nexus("mydata.nex")
#' pats <- quartet_site_pattern(aln, "human", "chimp", "gorilla", "orangutan")
#' pats$site_pattern[[1]]
#' }
quartet_site_pattern <- function(alignment, taxon1, taxon2, taxon3, taxon4,
                                 ...) {
  UseMethod("quartet_site_pattern")
}

#' @export
quartet_site_pattern.PhylotripletAlignment <- function(alignment, taxon1,
                                                       taxon2, taxon3, taxon4,
                                                       ...) {
  raw <- quartet_site_pattern_rcpp(alignment$.ptr, taxon1, taxon2, taxon3,
                                   taxon4)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    taxon4    = raw$taxon4,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- list(raw$site_pattern)
  df
}

#' @export
quartet_site_pattern.PhylotripletGene <- function(alignment, taxon1, taxon2,
                                                  taxon3, taxon4, ...) {
  raw <- gene_quartet_site_pattern_rcpp(alignment$.ptr, taxon1, taxon2,
                                        taxon3, taxon4)
  df <- data.frame(
    taxon1    = raw$taxon1,
    taxon2    = raw$taxon2,
    taxon3    = raw$taxon3,
    taxon4    = raw$taxon4,
    n_total   = raw$n_total,
    n_counted = raw$n_counted,
    stringsAsFactors = FALSE
  )
  df$site_pattern <- list(raw$site_pattern)
  df
}

#' Reshape a site-pattern vector into a 4^k array
#'
#' Converts a site-pattern vector of length 4^k (as stored in the
#' \code{site_pattern} column of \code{\link{triplet_site_patterns}} or the
#' \code{site_pattern} element of \code{\link{site_pattern}}) into a named
#' array of shape (4, 4, ..., 4) (k dimensions).  Element
#' \code{arr[b1, b2, ..., bk]} holds the fractional count of sites where
#' the first taxon shows base b1, the second b2, etc., with bases labelled
#' by the dimnames (A, C, G, T).
#'
#' @param vec Numeric vector of length 4^k.
#' @return A numeric array with k dimensions of size 4, with dimnames
#'   \code{c("A","C","G","T")} on each axis.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pats <- triplet_site_patterns(read_nexus("mydata.nex"))
#' arr  <- site_pattern_as_array(pats$site_pattern[[1]])
#' arr["A", "A", "A"]   # count for taxa showing A, A, A
#' }
site_pattern_as_array <- function(vec) {
  k <- log(length(vec), base = 4)
  if (abs(k - round(k)) > 1e-9 || k < 1)
    stop("vec length must be a power of 4 (got ", length(vec), ")")
  k <- as.integer(round(k))
  bases <- c("A", "C", "G", "T")
  # R fills arrays column-major, so array(vec, rep(4,k)) gives
  #   arr[a_k, ..., a_1] = vec[flat C-order index of (a_1,...,a_k)].
  # aperm with reversed permutation corrects axis order so that
  #   arr[a_1, a_2, ..., a_k] = vec[flat C-order index of (a_1,...,a_k)],
  # matching the C++ index b_0*4^(k-1) + b_1*4^(k-2) + ... + b_{k-1}.
  arr <- aperm(array(vec, dim = rep(4L, k)), rev(seq_len(k)))
  dimnames(arr) <- rep(list(bases), k)
  arr
}

#' Compute site-pattern counts for a k-tuple of taxa
#'
#' Accepts a vector of k taxon names (k >= 2), enabling pair, triplet,
#' quartet, or higher-order patterns.  Genes that do not contain all k taxa
#' are skipped with a warning.  Names may be given in any order.
#'
#' @param alignment A PhylotripletAlignment object.
#' @param taxa Character vector of k taxon names (any order, length >= 2).
#' @param ... Ignored.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{taxa}{Character vector of k taxon names in lexicographic order.}
#'     \item{n_total}{Total alignment columns across contributing genes.}
#'     \item{n_counted}{Columns with no gap in any of the k taxa.}
#'     \item{site_pattern}{Numeric vector of length 4^k with site-pattern
#'       counts.  Use \code{\link{site_pattern_as_array}} to reshape.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln  <- read_nexus("mydata.nex")
#' res  <- site_pattern(aln, c("human", "chimp", "gorilla"))
#' arr  <- site_pattern_as_array(res$site_pattern)
#' arr["A", "A", "A"]   # count for all three showing A
#' }
site_pattern <- function(alignment, taxa, ...) {
  UseMethod("site_pattern")
}

#' @export
site_pattern.PhylotripletAlignment <- function(alignment, taxa, ...) {
  taxa <- as.character(taxa)
  if (length(taxa) < 2L) stop("taxa must have length >= 2")
  raw <- site_pattern_rcpp(alignment$.ptr, taxa)
  list(
    taxa         = raw$taxa,
    n_total      = raw$n_total,
    n_counted    = raw$n_counted,
    site_pattern = raw$site_pattern
  )
}

#' @export
site_pattern.PhylotripletGene <- function(alignment, taxa, ...) {
  taxa <- as.character(taxa)
  if (length(taxa) < 2L) stop("taxa must have length >= 2")
  raw <- gene_site_pattern_rcpp(alignment$.ptr, taxa)
  list(
    taxa         = raw$taxa,
    n_total      = raw$n_total,
    n_counted    = raw$n_counted,
    site_pattern = raw$site_pattern
  )
}

#' Compute triplet topology weights
#'
#' Computes weights for the three possible rooted triplet topologies for a
#' specific taxon triplet.  The taxon names may be given in any order; they
#' are sorted lexicographically internally.
#'
#' The three returned weights correspond to topologies where a, b, c are the
#' taxon names in lexicographic order:
#' \describe{
#'   \item{w[1]}{(a,b)|c  -- a and b are sisters}
#'   \item{w[2]}{(a,c)|b  -- a and c are sisters}
#'   \item{w[3]}{(b,c)|a  -- b and c are sisters}
#' }
#'
#' @param alignment A PhylotripletAlignment or PhylotripletGene object.
#' @param taxon1,taxon2,taxon3 Character scalars naming the three taxa.
#' @param ... Ignored.
#'
#' @return A numeric vector of length 3 with weights for the three topologies.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln <- read_nexus("mydata.nex")
#' w <- compute_triplet_weights(aln, "human", "chimp", "gorilla")
#' }
compute_triplet_weights <- function(alignment, taxon1, taxon2, taxon3, ...) {
  UseMethod("compute_triplet_weights")
}

#' Compute JC69 triplet topology weights
#'
#' Computes likelihood-based weights for the three possible rooted triplet
#' topologies under the Jukes-Cantor (JC69) substitution model.
#'
#' @param alignment A \code{PhylotripletAlignment} or \code{PhylotripletGene}
#' @param taxon1,taxon2,taxon3 Taxon names (any order)
#' @param theta Population-scaled mutation rate parameter
#' @param ... Reserved for future use
#'
#' @return Named numeric vector of length 3: weights for topologies
#'   (ab|c, ac|b, bc|a) where a,b,c are sorted taxon names
#'
#' @export
compute_jc69_triplet_weights <- function(alignment, taxon1, taxon2, taxon3,
                                         theta, ...) {
  UseMethod("compute_jc69_triplet_weights")
}

#' @export
compute_jc69_triplet_weights.PhylotripletAlignment <- function(alignment,
                                                               taxon1, taxon2,
                                                               taxon3, theta,
                                                               ...) {
  compute_jc69_triplet_weights_rcpp(alignment$.ptr, taxon1, taxon2, taxon3,
                                    theta)
}

#' @export
compute_jc69_triplet_weights.PhylotripletGene <- function(alignment, taxon1,
                                                          taxon2, taxon3,
                                                          theta, ...) {
  gene_compute_jc69_triplet_weights_rcpp(alignment$.ptr, taxon1, taxon2,
                                         taxon3, theta)
}

#' GTR triplet topology weights
#'
#' Computes weights for the three rooted triplet topologies under the GTR
#' substitution model.
#'
#' @param alignment A \code{PhylotripletAlignment} or \code{PhylotripletGene}
#'   object.
#' @param taxon1,taxon2,taxon3 Taxon names (any order).
#' @param theta Population-scaled mutation rate parameter.
#' @param mu Mutation rate scaling parameter (default 1.0).
#' @param ... Ignored.
#'
#' @return Named numeric vector of length 3: weights for topologies
#'   (ab|c, ac|b, bc|a) where a,b,c are sorted taxon names
#'
#' @export
compute_gtr_triplet_weights <- function(alignment, taxon1, taxon2, taxon3,
                                        theta, mu = 1.0, ...) {
  UseMethod("compute_gtr_triplet_weights")
}

#' @export
compute_gtr_triplet_weights.PhylotripletAlignment <- function(alignment,
                                                              taxon1, taxon2,
                                                              taxon3, theta,
                                                              mu = 1.0, ...) {
  compute_gtr_triplet_weights_rcpp(alignment$.ptr, taxon1, taxon2, taxon3,
                                   theta, mu)
}

#' @export
compute_gtr_triplet_weights.PhylotripletGene <- function(alignment, taxon1,
                                                         taxon2, taxon3,
                                                         theta, mu = 1.0, ...) {
  gene_compute_gtr_triplet_weights_rcpp(alignment$.ptr, taxon1, taxon2,
                                        taxon3, theta, mu)
}

#' @export
compute_triplet_weights.PhylotripletAlignment <- function(alignment, taxon1,
                                                          taxon2, taxon3, ...) {
  compute_triplet_weights_rcpp(alignment$.ptr, taxon1, taxon2, taxon3)
}

#' @export
compute_triplet_weights.PhylotripletGene <- function(alignment, taxon1,
                                                     taxon2, taxon3, ...) {
  gene_compute_triplet_weights_rcpp(alignment$.ptr, taxon1, taxon2, taxon3)
}

#' Empirical nucleotide frequencies
#'
#' Returns empirical nucleotide frequencies (pi_A, pi_C, pi_G, pi_T) computed
#' from the alignment.  Gap columns are excluded; ambiguous bases contribute
#' fractionally (e.g. R counts as 0.5 A and 0.5 G).
#'
#' @param alignment A \code{PhylotripletAlignment} or \code{PhylotripletGene}
#'   object.
#' @param ... Ignored.
#'
#' @return A named numeric vector of length 4 with names
#'   \code{c("pi_A", "pi_C", "pi_G", "pi_T")} that sums to 1.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aln <- read_nexus("mydata.nex")
#' nucleotide_frequencies(aln)
#' }
nucleotide_frequencies <- function(alignment, ...) {
  UseMethod("nucleotide_frequencies")
}

#' @export
nucleotide_frequencies.PhylotripletAlignment <- function(alignment, ...) {
  nucleotide_frequencies_rcpp(alignment$.ptr)
}

#' @export
nucleotide_frequencies.PhylotripletGene <- function(alignment, ...) {
  gene_nucleotide_frequencies_rcpp(alignment$.ptr)
}