#pragma once
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "file_io.hpp"
#include "gen_tree.hpp"

namespace phylotriplet {

std::vector<WeightedTriplet> extract_weighted_triplets(
    const std::vector<WeightedTree>& weighted_trees) {
  std::vector<WeightedTriplet> weighted_triplets;
  for (const auto& wt : weighted_trees) {
    auto tree_triplets = wt.tree.get_triplets();
    for (const auto& t : tree_triplets) {
      weighted_triplets.emplace_back(t, wt.weight);
    }
  }
  return weighted_triplets;
}

std::string estimate_tree(const std::vector<WeightedTree>& weighted_trees) {
  auto weighted_triplets = extract_weighted_triplets(weighted_trees);
  PhyloTree tree = gen_tree(weighted_triplets);
  return tree.to_newick();
}

std::string estimate_tree_from_file(const std::string& filepath) {
  auto weighted_trees = read_weighted_trees_from_file(filepath);
  return estimate_tree(weighted_trees);
}

std::string estimate_tree_from_string(const std::string& content) {
  auto weighted_trees = read_weighted_trees_from_string(content);
  return estimate_tree(weighted_trees);
}

}  // namespace phylotriplet