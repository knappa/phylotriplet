#include <Rcpp.h>
#include <phylotriplet/triplet_estimator.hpp>
#include <phylotriplet/file_io.hpp>

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
		weighted_trees.emplace_back(
			as<std::string>(newicks[i]),
			weights[i]
		);
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

	return DataFrame::create(
		Named("newick") = newicks,
		Named("weight") = weights,
		Named("stringsAsFactors") = false
	);
}