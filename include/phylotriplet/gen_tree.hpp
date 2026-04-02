#pragma once
#include <Eigen/Dense>
#include <array>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "spectral_partition.hpp"
#include "tree.hpp"

namespace phylotriplet {

/**
 * A triplet with an associated weight.
 * triplet encodes a,b|c (a and b share a more recent common ancestor than c).
 */
struct WeightedTriplet {
  std::array<std::string, 3> triplet;
  double weight;

  WeightedTriplet(const std::array<std::string, 3>& trip)
      : triplet(trip), weight(1.0) {}
  WeightedTriplet(const std::array<std::string, 3>& trip, double w)
      : triplet(trip), weight(w) {}
};

namespace detail {

inline void gen_tree_recursive(const std::vector<WeightedTriplet>& triplets,
                               TreeNode* node) {
  // Base case: single triplet (a,b,c) encodes a,b|c.
  if (triplets.size() == 1) {
    const auto& t = triplets[0].triplet;
    node->add_child(t[2]);
    TreeNode* sub = node->add_child();
    sub->add_child(t[0]);
    sub->add_child(t[1]);
    return;
  }

  // Collect all taxa present in this set of triplets.
  std::set<std::string> taxa_set;
  for (const auto& wt : triplets) {
    taxa_set.insert(wt.triplet[0]);
    taxa_set.insert(wt.triplet[1]);
    taxa_set.insert(wt.triplet[2]);
  }
  std::vector<std::string> taxa(taxa_set.begin(), taxa_set.end());
  // std::set iteration is sorted, so taxa is already sorted.

  std::map<std::string, int> taxa_to_index;
  for (int i = 0; i < static_cast<int>(taxa.size()); ++i)
    taxa_to_index[taxa[i]] = i;

  // Build adjacency matrix: each triplet (a,b,c) contributes its weight
  // to the a–b edge (the cherry).
  const int n = static_cast<int>(taxa.size());
  Eigen::MatrixXd adj_matrix = Eigen::MatrixXd::Zero(n, n);
  for (const auto& wt : triplets) {
    const int ai = taxa_to_index[wt.triplet[0]];
    const int bi = taxa_to_index[wt.triplet[1]];
    adj_matrix(ai, bi) += wt.weight;
    adj_matrix(bi, ai) += wt.weight;
  }

  // Spectral Laplacian partition → two or more components.
  auto parts = spectral_laplacian_partition(adj_matrix, taxa);
  if (parts.size() < 2) {
    throw std::runtime_error(
        "gen_tree: spectral partition failed to split taxa; ");
  }

  for (const auto& comp : parts) {
    if (comp.empty()) {
      throw std::runtime_error(
          "gen_tree: spectral partition produced an empty component; this "
          "should not happen");
    } else if (comp.size() == 1) {
      // Single-taxon component: attach as a leaf.
      node->add_child(comp[0]);
    } else {
      // Multi-taxon component, we might have to recurse or emit a
      // multifurcation depending on the triplets over these taxa.

      // filter triplets to those over this component's taxa
      std::set<std::string> comp_set(comp.begin(), comp.end());
      std::vector<WeightedTriplet> comp_triplets;
      for (const auto& wt : triplets) {
        if (comp_set.count(wt.triplet[0]) && comp_set.count(wt.triplet[1]) &&
            comp_set.count(wt.triplet[2]))
          comp_triplets.push_back(wt);
      }

      TreeNode* subnode = node->add_child();
      if (comp_triplets.empty()) {
        // No triplets constrain this component: emit a multifurcation.
        for (const auto& m : comp)
          subnode->add_child(m);
      } else {
        gen_tree_recursive(comp_triplets, subnode);
      }
    }
  }
}

}  // namespace detail

/**
 * Generate a phylogenetic tree from a list of weighted rooted triplets.
 *
 * Each triplet is encoded as {a, b, c} meaning the topology a,b|c —
 * a and b share a more recent common ancestor than either does with c.
 * This is the same encoding used by PhyloTree::get_triplets().
 * Triplet weights influence the spectral partitioning step.
 *
 * The algorithm uses spectral Laplacian partitioning to recursively
 * bipartition the taxon set until each component is resolved into a
 * subtree.  When a component has taxa but no constraining triplets,
 * the taxa are attached as a multifurcation.
 *
 * @param triplets  Non-empty list of weighted triplets.
 * @return          A PhyloTree representing the inferred topology.
 * @throws std::invalid_argument if triplets is empty.
 * @throws std::runtime_error    if spectral partitioning cannot split.
 */
inline PhyloTree gen_tree(const std::vector<WeightedTriplet>& triplets) {
  if (triplets.empty())
    throw std::invalid_argument("gen_tree: triplets must be non-empty");

  PhyloTree tree;
  detail::gen_tree_recursive(triplets, tree.root.get());
  return tree;
}

}  // namespace phylotriplet
