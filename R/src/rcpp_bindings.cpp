#include <Rcpp.h>
#define PHYLOTRIPLET_RCPP_WARNINGS
#include <phylotriplet/file_io.hpp>
#include <phylotriplet/species_tree_jc69.hpp>
#include <phylotriplet/sequence_analysis.hpp>
#include <phylotriplet/sequence_io.hpp>
#include <phylotriplet/triplet_estimator.hpp>

using namespace Rcpp;
using namespace phylotriplet;

//' Estimate tree from weighted trees file (internal)
//'
//' @param filepath Path to file with weighted trees
//' @return Newick string
//' @keywords internal
// [[Rcpp::export]]
std::string estimate_tree_from_file_rcpp(std::string filepath) {
  return estimate_tree_from_file(filepath);
}

//' Estimate tree from weighted trees string (internal)
//'
//' @param content String containing weighted trees
//' @return Newick string
//' @keywords internal
// [[Rcpp::export]]
std::string estimate_tree_from_string_rcpp(std::string content) {
  return estimate_tree_from_string(content);
}

//' Estimate tree from data.frame (internal)
//'
//' @param trees_df Data frame with newick and weight columns
//' @return Newick string
//' @keywords internal
// [[Rcpp::export]]
std::string estimate_tree_rcpp(DataFrame trees_df) {
  CharacterVector newicks = trees_df["newick"];
  NumericVector weights = trees_df["weight"];

  std::vector<WeightedTree> weighted_trees;
  weighted_trees.reserve(trees_df.nrows());

  for (int i = 0; i < trees_df.nrows(); i++) {
    weighted_trees.emplace_back(as<std::string>(newicks[i]), weights[i]);
  }

  return estimate_tree(weighted_trees);
}

//' Read weighted trees from file (internal)
//'
//' @param filepath Path to weighted trees file
//' @return Data frame with newick and weight columns
//' @keywords internal
// [[Rcpp::export]]
DataFrame read_weighted_trees_rcpp(std::string filepath) {
  auto trees = read_weighted_trees_from_file(filepath);

  CharacterVector newicks(trees.size());
  NumericVector weights(trees.size());

  for (size_t i = 0; i < trees.size(); i++) {
    newicks[i] = trees[i].tree.to_newick();
    weights[i] = trees[i].weight;
  }

  return DataFrame::create(Named("newick") = newicks, Named("weight") = weights,
                           Named("stringsAsFactors") = false);
}

// -- Sequence alignment I/O
// ----------------------------------------------------

//' Read a NEXUS alignment file (internal)
//'
//' @param filepath Path to the NEXUS file
//' @return External pointer to a MultiGeneAlignment object
//' @keywords internal
// [[Rcpp::export]]
SEXP read_nexus_rcpp(std::string filepath) {
  MultiGeneAlignment* aln = new MultiGeneAlignment(read_nexus_file(filepath));
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln, true);
  return ptr;
}

//' Read one or more PHYLIP alignment files (internal)
//'
//' @param filepaths Character vector of file paths, one per gene
//' @param gene_names Optional character vector of gene names
//' @param interleaved Logical; use interleaved PHYLIP format
//' @return External pointer to a MultiGeneAlignment object
//' @keywords internal
// [[Rcpp::export]]
SEXP read_phylip_files_rcpp(std::vector<std::string> filepaths,
                            std::vector<std::string> gene_names,
                            bool interleaved) {
  MultiGeneAlignment* aln = new MultiGeneAlignment(
      read_phylip_files(filepaths, gene_names, interleaved));
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln, true);
  return ptr;
}

//' Get number of genes in a MultiGeneAlignment (internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @return Integer ' @keywords internal
// [[Rcpp::export]]
int alignment_num_genes_rcpp(SEXP aln_ptr) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  return static_cast<int>(ptr->num_genes());
}

//' Get sorted taxon names from a MultiGeneAlignment (internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @return Character vector of taxon names ' @keywords
// internal
// [[Rcpp::export]]
CharacterVector alignment_taxon_names_rcpp(SEXP aln_ptr) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto names = ptr->all_taxon_names();
  return CharacterVector(names.begin(), names.end());
}

//' Get gene names from a MultiGeneAlignment (internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @return Character vector of gene names in order '
//@keywords internal
// [[Rcpp::export]]
CharacterVector alignment_gene_names_rcpp(SEXP aln_ptr) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto names = ptr->gene_names();
  return CharacterVector(names.begin(), names.end());
}

//' Extract a single gene from a MultiGeneAlignment by index or name (internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @param idx 1-based integer index, or a character
// gene name ' @return External pointer to a GeneAlignment (copy) ' @keywords
// internal
// [[Rcpp::export]]
SEXP get_gene_rcpp(SEXP aln_ptr, SEXP idx) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  const auto& genes = ptr->genes;

  GeneAlignment* gene = nullptr;
  if (TYPEOF(idx) == INTSXP || TYPEOF(idx) == REALSXP) {
    int i = Rcpp::as<int>(idx) - 1;  // convert to 0-based
    if (i < 0 || static_cast<size_t>(i) >= genes.size())
      Rcpp::stop("gene index %d out of range (1..%d)", i + 1,
                 static_cast<int>(genes.size()));
    gene = new GeneAlignment(genes[static_cast<size_t>(i)]);
  } else if (TYPEOF(idx) == STRSXP) {
    std::string name = Rcpp::as<std::string>(idx);
    for (const auto& g : genes)
      if (g.gene_name == name) {
        gene = new GeneAlignment(g);
        break;
      }
    if (!gene) Rcpp::stop("gene '%s' not found", name.c_str());
  } else {
    Rcpp::stop("idx must be an integer index or a gene name string");
  }

  Rcpp::XPtr<GeneAlignment> gptr(gene, true);
  return gptr;
}

//' Get gene name from a GeneAlignment pointer (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Character scalar
//' @keywords internal
// [[Rcpp::export]]
std::string gene_name_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  return ptr->gene_name;
}

//' Get taxon names from a GeneAlignment pointer (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Character vector of taxon names
//' @keywords internal
// [[Rcpp::export]]
CharacterVector gene_taxon_names_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto names = ptr->taxon_names();
  return CharacterVector(names.begin(), names.end());
}

//' Check if a GeneAlignment has compressed site data (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Logical; TRUE if no compressed data is available
//' @keywords internal
// [[Rcpp::export]]
bool gene_compressed_empty_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  return ptr->compressed.empty();
}

//' Get the number of unique column patterns in a GeneAlignment (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Integer; 0 if not compressed
//' @keywords internal
// [[Rcpp::export]]
int gene_compressed_n_patterns_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  return static_cast<int>(ptr->compressed.n_patterns);
}

//' Get the alignment length (number of columns) for a GeneAlignment (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Integer
//' @keywords internal
// [[Rcpp::export]]
int gene_length_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  return static_cast<int>(ptr->length());
}

// -- Single-gene site-pattern functions
// ----------------------------------------

//' Compute all triplet site-pattern counts for one gene (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Named list with elements:
//'   taxon1, taxon2, taxon3 (character vectors),
//'   n_total, n_counted (numeric vectors),
//'   site_pattern (numeric matrix, one row per triplet, 64 columns)
//' @keywords internal
// [[Rcpp::export]]
List gene_triplet_site_patterns_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto patterns = triplet_site_patterns(*ptr);

  const int nrow = static_cast<int>(patterns.size());
  CharacterVector taxon1(nrow), taxon2(nrow), taxon3(nrow);
  NumericVector n_total_v(nrow), n_counted_v(nrow);
  NumericMatrix n_mat(nrow, 64);

  for (int i = 0; i < nrow; ++i) {
    const TuplePatterns& p = patterns[i];
    taxon1[i] = p.taxa[0];
    taxon2[i] = p.taxa[1];
    taxon3[i] = p.taxa[2];
    n_total_v[i] = p.n_total;
    n_counted_v[i] = p.n_counted;
    for (int k = 0; k < 64; ++k) n_mat(i, k) = p.site_pattern[k];
  }

  return List::create(Named("taxon1") = taxon1, Named("taxon2") = taxon2,
                      Named("taxon3") = taxon3, Named("n_total") = n_total_v,
                      Named("n_counted") = n_counted_v,
                      Named("site_pattern") = n_mat);
}

//' Compute site-pattern counts for a specific triplet in one gene (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @return Named list (same format as gene_triplet_site_patterns_rcpp, one row)
//' @keywords internal
// [[Rcpp::export]]
List gene_triplet_site_pattern_rcpp(SEXP gene_ptr, std::string taxon1,
                                    std::string taxon2, std::string taxon3) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  TuplePatterns p = triplet_site_pattern(*ptr, taxon1, taxon2, taxon3);
  NumericVector sp(64);
  for (int k = 0; k < 64; ++k) sp[k] = p.site_pattern[k];
  return List::create(Named("taxon1") = p.taxa[0], Named("taxon2") = p.taxa[1],
                      Named("taxon3") = p.taxa[2], Named("n_total") = p.n_total,
                      Named("n_counted") = p.n_counted,
                      Named("site_pattern") = sp);
}

//' Compute all quartet site-pattern counts for one gene (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Named list with elements:
//'   taxon1, taxon2, taxon3, taxon4 (character vectors),
//'   n_total, n_counted (numeric vectors),
//'   site_pattern (numeric matrix, one row per quartet, 256 columns)
//' @keywords internal
// [[Rcpp::export]]
List gene_quartet_site_patterns_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto patterns = quartet_site_patterns(*ptr);

  const int nrow = static_cast<int>(patterns.size());
  CharacterVector taxon1(nrow), taxon2(nrow), taxon3(nrow), taxon4(nrow);
  NumericVector n_total_v(nrow), n_counted_v(nrow);
  NumericMatrix n_mat(nrow, 256);

  for (int i = 0; i < nrow; ++i) {
    const TuplePatterns& p = patterns[i];
    taxon1[i] = p.taxa[0];
    taxon2[i] = p.taxa[1];
    taxon3[i] = p.taxa[2];
    taxon4[i] = p.taxa[3];
    n_total_v[i] = p.n_total;
    n_counted_v[i] = p.n_counted;
    for (int k = 0; k < 256; ++k) n_mat(i, k) = p.site_pattern[k];
  }

  return List::create(Named("taxon1") = taxon1, Named("taxon2") = taxon2,
                      Named("taxon3") = taxon3, Named("taxon4") = taxon4,
                      Named("n_total") = n_total_v,
                      Named("n_counted") = n_counted_v,
                      Named("site_pattern") = n_mat);
}

//' Compute site-pattern counts for a specific quartet in one gene (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxon1,taxon2,taxon3,taxon4 Taxon names (any order)
//' @return Named list (same format as gene_quartet_site_patterns_rcpp, one row)
//' @keywords internal
// [[Rcpp::export]]
List gene_quartet_site_pattern_rcpp(SEXP gene_ptr, std::string taxon1,
                                    std::string taxon2, std::string taxon3,
                                    std::string taxon4) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  TuplePatterns p = quartet_site_pattern(*ptr, taxon1, taxon2, taxon3, taxon4);
  NumericVector sp(256);
  for (int k = 0; k < 256; ++k) sp[k] = p.site_pattern[k];
  return List::create(Named("taxon1") = p.taxa[0], Named("taxon2") = p.taxa[1],
                      Named("taxon3") = p.taxa[2], Named("taxon4") = p.taxa[3],
                      Named("n_total") = p.n_total,
                      Named("n_counted") = p.n_counted,
                      Named("site_pattern") = sp);
}

//' Compute site-pattern counts for a k-tuple in one gene (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxa Character vector of taxon names (any order, length >= 2)
//' @return Named list: taxa, n_total, n_counted, site_pattern (4^k vector)
//' @keywords internal
// [[Rcpp::export]]
List gene_site_pattern_rcpp(SEXP gene_ptr, std::vector<std::string> taxa) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  TuplePatterns p = site_pattern(*ptr, taxa);
  const int sp_len = static_cast<int>(p.site_pattern.size());
  NumericVector sp(sp_len);
  for (int i = 0; i < sp_len; ++i) sp[i] = p.site_pattern[i];
  CharacterVector taxa_out(p.taxa.begin(), p.taxa.end());
  return List::create(Named("taxa") = taxa_out, Named("n_total") = p.n_total,
                      Named("n_counted") = p.n_counted,
                      Named("site_pattern") = sp);
}

//' Compute aggregate triplet site-pattern counts (internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @return Named list with elements: '   taxon1,
// taxon2, taxon3 (character vectors), '   n_total, n_counted (numeric vectors),
// ' site_pattern (numeric matrix, one row per triplet, 64 columns) ' @keywords
// internal
// [[Rcpp::export]]
List triplet_site_patterns_rcpp(SEXP aln_ptr) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto patterns = triplet_site_patterns(*ptr);

  const int nrow = static_cast<int>(patterns.size());
  CharacterVector taxon1(nrow), taxon2(nrow), taxon3(nrow);
  NumericVector n_total_v(nrow), n_counted_v(nrow);
  NumericMatrix n_mat(nrow, 64);

  for (int i = 0; i < nrow; ++i) {
    const TuplePatterns& p = patterns[i];
    taxon1[i] = p.taxa[0];
    taxon2[i] = p.taxa[1];
    taxon3[i] = p.taxa[2];
    n_total_v[i] = p.n_total;
    n_counted_v[i] = p.n_counted;
    for (int k = 0; k < 64; ++k) n_mat(i, k) = p.site_pattern[k];
  }

  return List::create(Named("taxon1") = taxon1, Named("taxon2") = taxon2,
                      Named("taxon3") = taxon3, Named("n_total") = n_total_v,
                      Named("n_counted") = n_counted_v,
                      Named("site_pattern") = n_mat);
}

//' Compute site-pattern counts for one specific taxon triplet (internal)
//'
//' Taxa names may be given in any order.  Genes that do not contain all
//' three taxa are skipped with a warning.
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @return Named list with elements:
//'   taxon1, taxon2, taxon3 (length-1 character vectors),
//'   n_total, n_counted (scalars),
//'   site_pattern (numeric vector of length 64)
//' @keywords internal
// [[Rcpp::export]]
List triplet_site_pattern_rcpp(SEXP aln_ptr, std::string taxon1,
                               std::string taxon2, std::string taxon3) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  TuplePatterns p = triplet_site_pattern(*ptr, taxon1, taxon2, taxon3);
  NumericVector sp(64);
  for (int k = 0; k < 64; ++k) sp[k] = p.site_pattern[k];
  return List::create(Named("taxon1") = p.taxa[0], Named("taxon2") = p.taxa[1],
                      Named("taxon3") = p.taxa[2], Named("n_total") = p.n_total,
                      Named("n_counted") = p.n_counted,
                      Named("site_pattern") = sp);
}

//' Compute aggregate quartet site-pattern counts (internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @return Named list with elements: '   taxon1,
// taxon2, taxon3, taxon4 (character vectors), '   n_total, n_counted (numeric
// vectors), '   site_pattern (numeric matrix, one row per quartet, 256 columns)
//' @keywords internal
// [[Rcpp::export]]
List quartet_site_patterns_rcpp(SEXP aln_ptr) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto patterns = quartet_site_patterns(*ptr);

  const int nrow = static_cast<int>(patterns.size());
  CharacterVector taxon1(nrow), taxon2(nrow), taxon3(nrow), taxon4(nrow);
  NumericVector n_total_v(nrow), n_counted_v(nrow);
  NumericMatrix n_mat(nrow, 256);

  for (int i = 0; i < nrow; ++i) {
    const TuplePatterns& p = patterns[i];
    taxon1[i] = p.taxa[0];
    taxon2[i] = p.taxa[1];
    taxon3[i] = p.taxa[2];
    taxon4[i] = p.taxa[3];
    n_total_v[i] = p.n_total;
    n_counted_v[i] = p.n_counted;
    for (int k = 0; k < 256; ++k) n_mat(i, k) = p.site_pattern[k];
  }

  return List::create(Named("taxon1") = taxon1, Named("taxon2") = taxon2,
                      Named("taxon3") = taxon3, Named("taxon4") = taxon4,
                      Named("n_total") = n_total_v,
                      Named("n_counted") = n_counted_v,
                      Named("site_pattern") = n_mat);
}

//' Compute site-pattern counts for one specific taxon quartet (internal)
//'
//' Taxa names may be given in any order.  Genes that do not contain all
//' four taxa are skipped with a warning.
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @param taxon1,taxon2,taxon3,taxon4 Taxon names (any
// order) ' @return Named list with elements: '   taxon1, taxon2, taxon3, taxon4
//(length-1 character vectors), '   n_total, n_counted (scalars), ' site_pattern
//(numeric vector of length 256) ' @keywords internal
// [[Rcpp::export]]
List quartet_site_pattern_rcpp(SEXP aln_ptr, std::string taxon1,
                               std::string taxon2, std::string taxon3,
                               std::string taxon4) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  TuplePatterns p = quartet_site_pattern(*ptr, taxon1, taxon2, taxon3, taxon4);
  NumericVector sp(256);
  for (int k = 0; k < 256; ++k) sp[k] = p.site_pattern[k];
  return List::create(Named("taxon1") = p.taxa[0], Named("taxon2") = p.taxa[1],
                      Named("taxon3") = p.taxa[2], Named("taxon4") = p.taxa[3],
                      Named("n_total") = p.n_total,
                      Named("n_counted") = p.n_counted,
                      Named("site_pattern") = sp);
}

//' Compute site-pattern counts for a k-tuple of taxa (internal)
//'
//' Taxa names may be given in any order (k >= 2).  Genes missing any taxon
//' are skipped with a warning.  Throws if a taxon is absent from every gene.
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @param taxa Character vector of taxon names (any
// order, length >= 2) ' @return Named list with elements: '   taxa (character
// vector of length k, sorted), '   n_total, n_counted (scalars), ' site_pattern
//(numeric vector of length 4^k) ' @keywords internal
// [[Rcpp::export]]
List site_pattern_rcpp(SEXP aln_ptr, std::vector<std::string> taxa) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  TuplePatterns p = site_pattern(*ptr, taxa);
  const int sp_len = static_cast<int>(p.site_pattern.size());
  NumericVector sp(sp_len);
  for (int i = 0; i < sp_len; ++i) sp[i] = p.site_pattern[i];
  CharacterVector taxa_out(p.taxa.begin(), p.taxa.end());
  return List::create(Named("taxa") = taxa_out, Named("n_total") = p.n_total,
                      Named("n_counted") = p.n_counted,
                      Named("site_pattern") = sp);
}

//' Compute JC69 triplet topology weights (multi-gene, internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_jc69_triplet_weights_rcpp(SEXP aln_ptr,
                                                std::string taxon1,
                                                std::string taxon2,
                                                std::string taxon3,
                                                double theta) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto w = compute_jc69_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute JC69 triplet topology weights (single gene, internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gene_compute_jc69_triplet_weights_rcpp(SEXP gene_ptr,
                                                     std::string taxon1,
                                                     std::string taxon2,
                                                     std::string taxon3,
                                                     double theta) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto w = compute_jc69_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute GTR triplet topology weights (multi-gene, internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp / read_phylip_files_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @param mu Mutation rate scaling parameter (default 1.0)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_gtr_triplet_weights_rcpp(SEXP aln_ptr,
                                               std::string taxon1,
                                               std::string taxon2,
                                               std::string taxon3,
                                               double theta,
                                               double mu = 1.0) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto w = compute_gtr_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta, mu);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute GTR triplet topology weights (single gene, internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @param mu Mutation rate scaling parameter (default 1.0)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gene_compute_gtr_triplet_weights_rcpp(SEXP gene_ptr,
                                                    std::string taxon1,
                                                    std::string taxon2,
                                                    std::string taxon3,
                                                    double theta,
                                                    double mu = 1.0) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto w = compute_gtr_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta, mu);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute triplet topology weights for a specific taxon triplet (multi-gene,
// internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
// read_phylip_files_rcpp ' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_triplet_weights_rcpp(SEXP aln_ptr, std::string taxon1,
                                           std::string taxon2,
                                           std::string taxon3) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto w = compute_triplet_weights(*ptr, taxon1, taxon2, taxon3);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute triplet topology weights for a specific taxon triplet (single gene,
// internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gene_compute_triplet_weights_rcpp(SEXP gene_ptr,
                                                std::string taxon1,
                                                std::string taxon2,
                                                std::string taxon3) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto w = compute_triplet_weights(*ptr, taxon1, taxon2, taxon3);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute empirical nucleotide frequencies for a multi-gene alignment
//' (internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp /
//'   read_phylip_files_rcpp
//' @return Named numeric vector of length 4: (pi_A, pi_C, pi_G, pi_T)
//' @keywords internal
// [[Rcpp::export]]
NumericVector nucleotide_frequencies_rcpp(SEXP aln_ptr) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto f = ptr->nucleotide_frequencies();
  NumericVector out =
      NumericVector::create(Named("pi_A") = f[0], Named("pi_C") = f[1],
                            Named("pi_G") = f[2], Named("pi_T") = f[3]);
  return out;
}

//' Compute empirical nucleotide frequencies for a single gene (internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @return Named numeric vector of length 4: (pi_A, pi_C, pi_G, pi_T)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gene_nucleotide_frequencies_rcpp(SEXP gene_ptr) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto f = ptr->nucleotide_frequencies();
  NumericVector out =
      NumericVector::create(Named("pi_A") = f[0], Named("pi_C") = f[1],
                            Named("pi_G") = f[2], Named("pi_T") = f[3]);
  return out;
}

//' Compute JC69 species tree triplet topology weights (multi-gene, internal)
//'
//' @param aln_ptr External pointer from read_nexus_rcpp / read_phylip_files_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @param mu JC69 branch-length scaling (default 4/3)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_jc69_species_triplet_weights_rcpp(SEXP aln_ptr,
                                                        std::string taxon1,
                                                        std::string taxon2,
                                                        std::string taxon3,
                                                        double theta,
                                                        double mu = 4.0/3.0) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto w = compute_jc69_species_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta, mu);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute K2P species tree triplet topology weights (multi-gene, internal)
//'
//' @param aln_ptr External pointer from read_alignment_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @param mu K2P total substitution rate alpha+2*beta (default 4/3)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_k2p_species_triplet_weights_rcpp(SEXP aln_ptr,
                                                       std::string taxon1,
                                                       std::string taxon2,
                                                       std::string taxon3,
                                                       double theta,
                                                       double mu = 4.0/3.0) {
  Rcpp::XPtr<MultiGeneAlignment> ptr(aln_ptr);
  auto w = compute_k2p_species_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta, mu);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute JC69 species tree triplet topology weights (single gene, internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @param mu JC69 branch-length scaling (default 4/3)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gene_compute_jc69_species_triplet_weights_rcpp(SEXP gene_ptr,
                                                             std::string taxon1,
                                                             std::string taxon2,
                                                             std::string taxon3,
                                                             double theta,
                                                             double mu = 4.0/3.0) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto w = compute_jc69_species_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta, mu);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}

//' Compute K2P species tree triplet topology weights (single gene, internal)
//'
//' @param gene_ptr External pointer from get_gene_rcpp
//' @param taxon1,taxon2,taxon3 Taxon names (any order)
//' @param theta Population-scaled mutation rate parameter
//' @param mu K2P total substitution rate alpha+2*beta (default 4/3)
//' @return Numeric vector of length 3: weights for the three rooted triplet
//'   topologies (ab|c, ac|b, bc|a where a,b,c are the sorted taxon names)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gene_compute_k2p_species_triplet_weights_rcpp(SEXP gene_ptr,
                                                            std::string taxon1,
                                                            std::string taxon2,
                                                            std::string taxon3,
                                                            double theta,
                                                            double mu = 4.0/3.0) {
  Rcpp::XPtr<GeneAlignment> ptr(gene_ptr);
  auto w = compute_k2p_species_triplet_weights(*ptr, taxon1, taxon2, taxon3, theta, mu);
  NumericVector out(3);
  out[0] = w[0];
  out[1] = w[1];
  out[2] = w[2];
  return out;
}
