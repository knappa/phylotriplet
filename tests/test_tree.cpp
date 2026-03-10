#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "phylotriplet/tree.hpp"

using namespace phylotriplet;

// --- helpers -----------------------------------------------------------------

static std::set<std::string> leaf_set(const PhyloTree& t) {
  auto v = t.get_leaf_names();
  return std::set<std::string>(v.begin(), v.end());
}

static std::set<std::array<std::string, 3>> triplet_set(const PhyloTree& t) {
  auto v = t.get_triplets();
  return std::set<std::array<std::string, 3>>(v.begin(), v.end());
}

// --- Newick parse / emit
// ------------------------------------------------------

TEST_CASE("Newick round-trip: simple cherry") {
  PhyloTree t = PhyloTree::from_newick("(A,B);");
  REQUIRE(leaf_set(t) == std::set<std::string>{"A", "B"});
  // Both orderings are valid Newick
  std::string nwk = t.to_newick();
  REQUIRE((nwk == "(A,B);" || nwk == "(B,A);"));
}

TEST_CASE("Newick round-trip: three taxa") {
  std::string input = "((A,B),C);";
  PhyloTree t = PhyloTree::from_newick(input);
  REQUIRE(leaf_set(t) == std::set<std::string>{"A", "B", "C"});
  REQUIRE(t.to_newick() == input);
}

TEST_CASE("Newick round-trip: branch lengths preserved") {
  PhyloTree t = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);");
  REQUIRE(leaf_set(t) == std::set<std::string>{"A", "B", "C"});
  // Check branch lengths survive the round-trip
  std::string nwk = t.to_newick();
  REQUIRE(nwk.find("0.1") != std::string::npos);
  REQUIRE(nwk.find("0.4") != std::string::npos);
}

TEST_CASE("Newick round-trip: four-taxon symmetric tree") {
  std::string input = "((A,B),(C,D));";
  PhyloTree t = PhyloTree::from_newick(input);
  REQUIRE(leaf_set(t) == std::set<std::string>{"A", "B", "C", "D"});
  REQUIRE(t.to_newick() == input);
}

TEST_CASE("Newick: single leaf") {
  PhyloTree t = PhyloTree::from_newick("A;");
  REQUIRE(leaf_set(t) == std::set<std::string>{"A"});
  REQUIRE(t.to_newick() == "A;");
}

// --- get_leaf_names
// -----------------------------------------------------------

TEST_CASE("get_leaf_names: five taxa") {
  PhyloTree t = PhyloTree::from_newick("(((A,B),C),(D,E));");
  REQUIRE(leaf_set(t) == std::set<std::string>{"A", "B", "C", "D", "E"});
}

// --- get_triplets
// -------------------------------------------------------------

TEST_CASE("get_triplets: three-taxon tree has exactly one triplet") {
  // ((A,B),C) -> only triplet is {A,B,C} (A,B form the cherry)
  PhyloTree t = PhyloTree::from_newick("((A,B),C);");
  auto trips = t.get_triplets();
  REQUIRE(trips.size() == 1);
  std::array<std::string, 3> trip = trips[0];
  // cherry = {A,B}, outgroup = C
  REQUIRE(((trip[0] == "A" || trip[0] == "B") &&
           (trip[1] == "A" || trip[1] == "B") && trip[0] != trip[1] &&
           trip[2] == "C"));
}

TEST_CASE("get_triplets: four-taxon symmetric tree") {
  // ((A,B),(C,D))
  // triplets: {A,B | C}, {A,B | D}, {C,D | A}, {C,D | B}
  PhyloTree t = PhyloTree::from_newick("((A,B),(C,D));");
  auto trips = triplet_set(t);
  REQUIRE(trips.size() == 4);
  REQUIRE(trips.count({"A", "B", "C"}));
  REQUIRE(trips.count({"A", "B", "D"}));
  REQUIRE(trips.count({"C", "D", "A"}));
  REQUIRE(trips.count({"C", "D", "B"}));
}

TEST_CASE("get_triplets: caterpillar (((A,B),C),D)") {
  // Cherries at each level:
  //   node(A,B): cherry {A,B}, outgroups C, D -> {A,B,C} {A,B,D}
  //   node((A,B),C): cherry {A,B,C subtree}, outgroup D
  //     pairs from {A,B} x {D} and {B} x {D} and {A} x {D}:
  //     {A,C,D} {B,C,D}
  // Total: 4 triplets
  PhyloTree t = PhyloTree::from_newick("(((A,B),C),D);");
  auto trips = triplet_set(t);
  REQUIRE(trips.size() == 4);
  REQUIRE(trips.count({"A", "B", "C"}));
  REQUIRE(trips.count({"A", "B", "D"}));
  REQUIRE(trips.count({"A", "C", "D"}));
  REQUIRE(trips.count({"B", "C", "D"}));
}

// --- Robinson-Foulds
// ----------------------------------------------------------

TEST_CASE("robinson_foulds: identical trees") {
  PhyloTree t1 = PhyloTree::from_newick("((A,B),(C,D));");
  PhyloTree t2 = PhyloTree::from_newick("((A,B),(C,D));");
  REQUIRE(t1.robinson_foulds(t2) == 0);
}

TEST_CASE("robinson_foulds: three-taxon trees are always identical topology") {
  PhyloTree t1 = PhyloTree::from_newick("((A,B),C);");
  PhyloTree t2 = PhyloTree::from_newick("((A,B),C);");
  REQUIRE(t1.robinson_foulds(t2) == 0);
}

TEST_CASE("robinson_foulds: maximally different 4-taxon trees") {
  // Each tree has one non-trivial split. Neither split is shared.
  // RF = 1 (only in t1) + 1 (only in t2) = 2.
  PhyloTree t1 = PhyloTree::from_newick("((A,B),(C,D));");
  PhyloTree t2 = PhyloTree::from_newick("((A,C),(B,D));");
  REQUIRE(t1.robinson_foulds(t2) == 2);
}

TEST_CASE("robinson_foulds: symmetric") {
  PhyloTree t1 = PhyloTree::from_newick("((A,B),(C,D));");
  PhyloTree t2 = PhyloTree::from_newick("((A,C),(B,D));");
  REQUIRE(t1.robinson_foulds(t2) == t2.robinson_foulds(t1));
}

// --- prune -------------------------------------------------------------------

TEST_CASE("prune: retains only requested taxa") {
  PhyloTree t = PhyloTree::from_newick("((A,B),(C,D));");
  t.prune({"A", "C", "D"});
  REQUIRE(leaf_set(t) == std::set<std::string>{"A", "C", "D"});
}

TEST_CASE("prune: to three taxa yields correct topology") {
  // ((A,B),(C,D)) pruned to {A,C,D} -> ((C,D),A) or equivalent
  PhyloTree t = PhyloTree::from_newick("((A,B),(C,D));");
  t.prune({"A", "C", "D"});
  // The C,D cherry must still exist
  auto trips = triplet_set(t);
  REQUIRE(trips.size() == 1);
  REQUIRE(trips.count({"C", "D", "A"}));
}

TEST_CASE("prune: to two taxa") {
  PhyloTree t = PhyloTree::from_newick("(((A,B),C),(D,E));");
  t.prune({"B", "E"});
  REQUIRE(leaf_set(t) == std::set<std::string>{"B", "E"});
}

// --- LCA ---------------------------------------------------------------------

TEST_CASE("get_common_ancestor: cherry returns internal node") {
  PhyloTree t = PhyloTree::from_newick("((A,B),C);");
  TreeNode* lca = t.root->get_common_ancestor({"A", "B"});
  REQUIRE(lca != nullptr);
  REQUIRE(!lca->is_leaf());
  auto names = lca->get_leaf_names();
  REQUIRE(std::set<std::string>(names.begin(), names.end()) ==
          std::set<std::string>{"A", "B"});
}

TEST_CASE("get_common_ancestor: all taxa returns root") {
  PhyloTree t = PhyloTree::from_newick("((A,B),C);");
  TreeNode* lca = t.root->get_common_ancestor({"A", "B", "C"});
  REQUIRE(lca == t.root.get());
}

TEST_CASE("get_common_ancestor: single taxon returns that leaf") {
  PhyloTree t = PhyloTree::from_newick("((A,B),C);");
  TreeNode* lca = t.root->get_common_ancestor({"A"});
  REQUIRE(lca != nullptr);
  REQUIRE(lca->is_leaf());
  REQUIRE(lca->name == "A");
}

TEST_CASE("get_common_ancestor: missing taxon throws") {
  PhyloTree t = PhyloTree::from_newick("((A,B),C);");
  REQUIRE_THROWS_AS(t.root->get_common_ancestor({"X"}), std::runtime_error);
}

// --- traversal ---------------------------------------------------------------

TEST_CASE("traverse postorder: leaves come before their parent") {
  PhyloTree t = PhyloTree::from_newick("((A,B),C);");
  auto nodes = t.root->traverse("postorder");
  // In postorder, both A and B must appear before the (A,B) internal node,
  // and (A,B) internal must appear before the root.
  std::map<TreeNode*, size_t> pos;
  for (size_t i = 0; i < nodes.size(); ++i) pos[nodes[i]] = i;

  for (auto* n : nodes) {
    for (const auto& child : n->children) {
      REQUIRE(pos[child.get()] < pos[n]);
    }
  }
}

TEST_CASE("traverse preorder: parent comes before its children") {
  PhyloTree t = PhyloTree::from_newick("((A,B),C);");
  auto nodes = t.root->traverse("preorder");
  std::map<TreeNode*, size_t> pos;
  for (size_t i = 0; i < nodes.size(); ++i) pos[nodes[i]] = i;

  for (auto* n : nodes) {
    for (const auto& child : n->children) {
      REQUIRE(pos[n] < pos[child.get()]);
    }
  }
}

TEST_CASE("traverse levelorder: depth is non-decreasing") {
  PhyloTree t = PhyloTree::from_newick("(((A,B),C),(D,E));");
  auto nodes = t.root->traverse("levelorder");

  // Compute depth of each node.
  std::map<TreeNode*, int> depth;
  depth[t.root.get()] = 0;
  for (auto* n : nodes) {
    for (const auto& child : n->children) {
      depth[child.get()] = depth[n] + 1;
    }
  }

  for (size_t i = 1; i < nodes.size(); ++i) {
    REQUIRE(depth[nodes[i - 1]] <= depth[nodes[i]]);
  }
}

TEST_CASE("traverse: total node count matches tree size") {
  // ((A,B),(C,D)) has 4 leaves + 3 internal nodes = 7 nodes
  PhyloTree t = PhyloTree::from_newick("((A,B),(C,D));");
  REQUIRE(t.root->traverse("postorder").size() == 7);
  REQUIRE(t.root->traverse("preorder").size() == 7);
  REQUIRE(t.root->traverse("levelorder").size() == 7);
}
