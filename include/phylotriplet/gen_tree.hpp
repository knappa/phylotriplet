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

  WeightedTriplet(const std::array<std::string, 3>& trip) : triplet(trip), weight(1.0) {}
  WeightedTriplet(const std::array<std::string, 3>& trip, double w) : triplet(trip), weight(w) {}
};

namespace detail {

/**
 * Resolve `taxa` into a subtree hanging below `node`, using `triplets`.
 *
 * `taxa` is the authoritative taxon set for this node: every name in it ends
 * up as a leaf of the emitted subtree.  `triplets` only supplies topological
 * constraints, and may mention a strict subset of `taxa` (a taxon whose
 * triplets all straddle the current partition boundary contributes no
 * constraint here, but must still be placed).
 */
inline void gen_tree_recursive(const std::vector<WeightedTriplet>& triplets,
                               const std::vector<std::string>& taxa,
                               TreeNode* node) {
  // One or two taxa admit no further resolution: attach them directly.
  if (taxa.size() <= 2) {
    for (const auto& name : taxa)
      node->add_child(name);
    return;
  }

  // No constraints over these taxa: emit a multifurcation.
  if (triplets.empty()) {
    for (const auto& name : taxa)
      node->add_child(name);
    return;
  }

  // Base case: a single triplet (a,b,c) over exactly these taxa encodes a,b|c.
  if (triplets.size() == 1 && taxa.size() == 3) {
    const auto& t = triplets[0].triplet;
    node->add_child(t[2]);
    TreeNode* sub = node->add_child();
    sub->add_child(t[0]);
    sub->add_child(t[1]);
    return;
  }

  std::map<std::string, int> taxa_to_index;
  for (int i = 0; i < static_cast<int>(taxa.size()); ++i)
    taxa_to_index[taxa[i]] = i;

  // Build adjacency matrix: each triplet (a,b,c) contributes its weight
  // to the a-b edge (the cherry).  A taxon mentioned by no triplet stays an
  // isolated vertex, which the partition reports as its own component.
  const int n = static_cast<int>(taxa.size());
  Eigen::MatrixXd adj_matrix = Eigen::MatrixXd::Zero(n, n);
  for (const auto& wt : triplets) {
    const int ai = taxa_to_index[wt.triplet[0]];
    const int bi = taxa_to_index[wt.triplet[1]];
    adj_matrix(ai, bi) += wt.weight;
    adj_matrix(bi, ai) += wt.weight;
  }

  // Spectral Laplacian partition -> two or more components.
  auto parts = spectral_laplacian_partition(adj_matrix, taxa);
  if (parts.size() < 2) {
    throw std::runtime_error("gen_tree: spectral partition failed to split taxa; ");
  }

  for (const auto& comp : parts) {
    if (comp.empty()) {
      throw std::runtime_error(
          "gen_tree: spectral partition produced an empty component; this "
          "should not happen");
    } else if (comp.size() == 1) {
      // Single-taxon component: attach as a leaf.
      node->add_child(comp[0]);
    } else if (comp.size() == taxa.size()) {
      // A component holding every taxon is not a split, and recursing on it
      // would not terminate.  The partition cannot report one, so reaching
      // here is a bug in the partitioning code rather than bad input.
      throw std::runtime_error(
          "gen_tree: spectral partition returned a component containing every "
          "taxon; this should not happen");
    } else {
      // Multi-taxon component: recurse on the triplets wholly inside it.
      std::set<std::string> comp_set(comp.begin(), comp.end());
      std::vector<WeightedTriplet> comp_triplets;
      for (const auto& wt : triplets) {
        if (comp_set.count(wt.triplet[0]) && comp_set.count(wt.triplet[1]) &&
            comp_set.count(wt.triplet[2]))
          comp_triplets.push_back(wt);
      }

      TreeNode* subnode = node->add_child();
      gen_tree_recursive(comp_triplets, comp, subnode);
    }
  }
}

}  // namespace detail

/**
 * Generate a phylogenetic tree from a list of weighted rooted triplets.
 *
 * Each triplet is encoded as {a, b, c} meaning the topology a,b|c:
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

  std::set<std::string> taxa_set;
  for (const auto& wt : triplets) {
    taxa_set.insert(wt.triplet[0]);
    taxa_set.insert(wt.triplet[1]);
    taxa_set.insert(wt.triplet[2]);
  }
  std::vector<std::string> taxa(taxa_set.begin(), taxa_set.end());

  PhyloTree tree;
  detail::gen_tree_recursive(triplets, taxa, tree.root.get());
  return tree;
}

}  // namespace phylotriplet
