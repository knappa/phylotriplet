#pragma once
#include <string>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <queue>
#include <cctype>

namespace phylotriplet {

/**
 * A node in a phylogenetic tree.
 *
 * Each node is owned by its parent via shared_ptr; the parent pointer is
 * non-owning (raw) to avoid cycles.  Leaf nodes are identified by having no
 * children.
 */
class TreeNode {
public:
	std::string name;       // leaf name, or internal node label
	double      dist = 0.0; // branch length to parent

	TreeNode* parent = nullptr;
	std::vector<std::shared_ptr<TreeNode>> children;

	TreeNode() = default;
	explicit TreeNode(std::string n, double d = 0.0)
		: name(std::move(n)), dist(d) {}

	bool is_leaf() const { return children.empty(); }
	bool is_root() const { return parent == nullptr; }

	// -----------------------------------------------------------------------
	// Children management
	// -----------------------------------------------------------------------

	/**
	 * Attach an existing node as a child.  Optionally override its dist.
	 * Returns a raw pointer to the newly-added child (owned by this node).
	 */
	TreeNode* add_child(std::shared_ptr<TreeNode> child, double d = -1.0) {
		if (d >= 0.0) child->dist = d;
		child->parent = this;
		children.push_back(std::move(child));
		return children.back().get();
	}

	/**
	 * Create a new leaf child with the given name and branch length.
	 * Returns a raw pointer to the newly-added child.
	 */
	TreeNode* add_child(const std::string& n = "", double d = 0.0) {
		return add_child(std::make_shared<TreeNode>(n, d));
	}

	/**
	 * Remove this node from its parent and return the shared_ptr that held it.
	 * Returns nullptr if this node has no parent (i.e. is the root).
	 */
	std::shared_ptr<TreeNode> detach() {
		if (!parent) return nullptr;
		auto& siblings = parent->children;
		for (auto it = siblings.begin(); it != siblings.end(); ++it) {
			if (it->get() == this) {
				auto self = *it;
				siblings.erase(it);
				parent = nullptr;
				return self;
			}
		}
		return nullptr;
	}

	// -----------------------------------------------------------------------
	// Traversal
	// -----------------------------------------------------------------------

	/**
	 * Return all nodes reachable from this node in the requested order.
	 * @param strategy  "postorder" (default), "preorder", or "levelorder"
	 */
	std::vector<TreeNode*> traverse(const std::string& strategy = "postorder") const {
		std::vector<TreeNode*> result;
		if (strategy == "preorder") {
			traverse_preorder(result);
		} else if (strategy == "levelorder") {
			traverse_levelorder(result);
		} else {
			traverse_postorder(result);
		}
		return result;
	}

	/** Return all leaf nodes under this node. */
	std::vector<TreeNode*> get_leaves() const {
		std::vector<TreeNode*> leaves;
		for (auto* node : traverse()) {
			if (node->is_leaf()) leaves.push_back(node);
		}
		return leaves;
	}

	/** Return the names of all leaf nodes under this node, in postorder. */
	std::vector<std::string> get_leaf_names() const {
		std::vector<std::string> names;
		for (auto* leaf : get_leaves()) names.push_back(leaf->name);
		return names;
	}

	// -----------------------------------------------------------------------
	// Topology queries
	// -----------------------------------------------------------------------

	/**
	 * Find the lowest common ancestor of the leaves named in `taxa`.
	 * Throws std::runtime_error if any name is not present as a leaf.
	 */
	TreeNode* get_common_ancestor(const std::vector<std::string>& taxa) const {
		if (taxa.empty()) return nullptr;

		std::map<std::string, TreeNode*> leaf_map;
		for (auto* node : traverse()) {
			if (node->is_leaf()) leaf_map[node->name] = node;
		}

		std::vector<TreeNode*> nodes;
		for (const auto& t : taxa) {
			auto it = leaf_map.find(t);
			if (it == leaf_map.end())
				throw std::runtime_error("Taxon not found in tree: " + t);
			nodes.push_back(it->second);
		}
		return lca_of(nodes);
	}

	// -----------------------------------------------------------------------
	// Newick serialisation
	// -----------------------------------------------------------------------

	/** Serialise this subtree as a Newick string (includes trailing semicolon). */
	std::string to_newick() const {
		std::string s;
		newick_recursive(s);
		s += ';';
		return s;
	}

	// -----------------------------------------------------------------------
	// Pruning
	// -----------------------------------------------------------------------

	/**
	 * Retain only the taxa listed in `taxa_to_keep` and contract any internal
	 * nodes that become degree-2 as a result (their branch lengths are summed).
	 * Does nothing to the root even if it ends up with one child.
	 */
	void prune(const std::vector<std::string>& taxa_to_keep) {
		std::set<std::string> keep(taxa_to_keep.begin(), taxa_to_keep.end());

		// Remove unwanted leaves, one at a time (safe with raw-pointer traversal).
		bool changed = true;
		while (changed) {
			changed = false;
			for (auto* node : traverse()) {
				if (node->is_leaf() && !node->name.empty()
						&& !keep.count(node->name) && node->parent) {
					auto& sib = node->parent->children;
					sib.erase(std::remove_if(sib.begin(), sib.end(),
						[node](const std::shared_ptr<TreeNode>& p) {
							return p.get() == node;
						}), sib.end());
					changed = true;
					break; // restart with fresh traversal
				}
			}
		}

		// Contract non-root internal nodes that now have exactly one child.
		changed = true;
		while (changed) {
			changed = false;
			for (auto* node : traverse()) {
				if (!node->is_leaf() && !node->is_root()
						&& node->children.size() == 1) {
					auto child = node->children[0];
					child->dist += node->dist;
					TreeNode* par = node->parent;
					// Replace node with child in par's children list.
					for (auto& slot : par->children) {
						if (slot.get() == node) {
							child->parent = par;
							slot = child;
							break;
						}
					}
					changed = true;
					break;
				}
			}
		}
	}

private:
	void traverse_preorder(std::vector<TreeNode*>& out) const {
		out.push_back(const_cast<TreeNode*>(this));
		for (const auto& c : children) c->traverse_preorder(out);
	}

	void traverse_postorder(std::vector<TreeNode*>& out) const {
		for (const auto& c : children) c->traverse_postorder(out);
		out.push_back(const_cast<TreeNode*>(this));
	}

	void traverse_levelorder(std::vector<TreeNode*>& out) const {
		std::queue<const TreeNode*> q;
		q.push(this);
		while (!q.empty()) {
			const TreeNode* n = q.front(); q.pop();
			out.push_back(const_cast<TreeNode*>(n));
			for (const auto& c : n->children) q.push(c.get());
		}
	}

	/** Lowest common ancestor of a non-empty set of nodes. */
	static TreeNode* lca_of(const std::vector<TreeNode*>& nodes) {
		if (nodes.size() == 1) return nodes[0];

		// Ancestor sets as depth-indexed: collect ancestors of first node.
		std::set<TreeNode*> ancestors;
		for (TreeNode* p = nodes[0]; p; p = p->parent) ancestors.insert(p);

		// Intersect with ancestors of each subsequent node.
		for (size_t i = 1; i < nodes.size(); ++i) {
			std::set<TreeNode*> cur_anc;
			for (TreeNode* p = nodes[i]; p; p = p->parent) cur_anc.insert(p);
			std::set<TreeNode*> inter;
			for (TreeNode* a : ancestors)
				if (cur_anc.count(a)) inter.insert(a);
			ancestors = std::move(inter);
		}

		// The LCA is the deepest (maximum depth) node in the intersection.
		TreeNode* lca = nullptr;
		int max_depth = -1;
		for (TreeNode* a : ancestors) {
			int d = 0;
			for (TreeNode* p = a->parent; p; p = p->parent) ++d;
			if (d > max_depth) { max_depth = d; lca = a; }
		}
		return lca;
	}

	/**
	 * Quote a Newick name if it contains characters that require quoting
	 * (spaces, parentheses, comma, colon, semicolon, or single quote).
	 * Embedded single quotes are escaped by doubling them per the Newick spec.
	 */
	static std::string quote_newick_name(const std::string& n) {
		bool needs_quotes = false;
		for (char c : n) {
			if (std::isspace(static_cast<unsigned char>(c)) ||
					c == '(' || c == ')' || c == ',' ||
					c == ':' || c == ';' || c == '\'') {
				needs_quotes = true;
				break;
			}
		}
		if (!needs_quotes) return n;
		std::string result = "'";
		for (char c : n) {
			if (c == '\'') result += "''";
			else result += c;
		}
		result += "'";
		return result;
	}

	void newick_recursive(std::string& s) const {
		if (!children.empty()) {
			s += '(';
			for (size_t i = 0; i < children.size(); ++i) {
				if (i) s += ',';
				children[i]->newick_recursive(s);
			}
			s += ')';
			if (!name.empty()) s += quote_newick_name(name);
		} else {
			s += quote_newick_name(name);
		}
		if (dist != 0.0) {
			std::ostringstream oss;
			oss << ':' << dist;
			s += oss.str();
		}
	}
};

// ---------------------------------------------------------------------------
// Newick parser (recursive descent)
// ---------------------------------------------------------------------------
namespace detail {

inline void skip_ws(const std::string& s, size_t& pos) {
	while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
}

// Parse a bare or single-, or double-quoted name token.
inline std::string parse_name(const std::string& s, size_t& pos) {
	std::string name;
	if (pos < s.size() && s[pos] == '\'') {
		++pos;
		while (pos < s.size() && s[pos] != '\'') name += s[pos++];
		if (pos < s.size()) ++pos;
		return name;
	} else if (pos < s.size() && s[pos] == '\"') {
		++pos;
		while (pos < s.size() && s[pos] != '\"') name += s[pos++];
		if (pos < s.size()) ++pos;
		return name;
	}
	while (pos < s.size()
			&& s[pos] != '(' && s[pos] != ')' && s[pos] != ','
			&& s[pos] != ':' && s[pos] != ';'
			&& !std::isspace(static_cast<unsigned char>(s[pos]))) {
		name += s[pos++];
	}
	return name;
}

inline double parse_double(const std::string& s, size_t& pos) {
	std::string num;
	if (pos < s.size() && (s[pos] == '-' || s[pos] == '+')) num += s[pos++];
	while (pos < s.size()
			&& (std::isdigit(static_cast<unsigned char>(s[pos]))
				|| s[pos] == '.' || s[pos] == 'e' || s[pos] == 'E'
				|| s[pos] == '-' || s[pos] == '+')) {
		num += s[pos++];
	}
	if (num.empty() || num == "-" || num == "+")
		throw std::invalid_argument(
			"Newick parse error: expected a numeric branch length at position "
			+ std::to_string(pos));
	return std::stod(num);
}

// Forward declaration.
inline std::shared_ptr<TreeNode> parse_node(const std::string& s, size_t& pos);

inline std::shared_ptr<TreeNode> parse_node(const std::string& s, size_t& pos) {
	skip_ws(s, pos);
	auto node = std::make_shared<TreeNode>();

	if (pos < s.size() && s[pos] == '(') {
		++pos; // consume '('
		while (pos < s.size() && s[pos] != ')') {
			skip_ws(s, pos);
			auto child = parse_node(s, pos);
			child->parent = node.get();
			node->children.push_back(std::move(child));
			skip_ws(s, pos);
			if (pos < s.size() && s[pos] == ',') ++pos;
		}
		if (pos < s.size()) ++pos; // consume ')'

		// Optional label (internal node name; ignored if numeric).
		skip_ws(s, pos);
		if (pos < s.size() && s[pos] != ':' && s[pos] != ','
				&& s[pos] != ')' && s[pos] != ';') {
			node->name = parse_name(s, pos);
		}
	} else {
		node->name = parse_name(s, pos);
	}

	// Optional branch length.
	skip_ws(s, pos);
	if (pos < s.size() && s[pos] == ':') {
		++pos;
		node->dist = parse_double(s, pos);
	}
	return node;
}

} // namespace detail

// ---------------------------------------------------------------------------
// PhyloTree — thin wrapper around a root TreeNode
// ---------------------------------------------------------------------------

/**
 * A rooted phylogenetic tree.
 *
 * Newick round-trip:
 * @code
 *   auto t = PhyloTree::from_newick("((human,chimp),gorilla);");
 *   std::string nwk = t.to_newick();
 * @endcode
 *
 * Robinson-Foulds distance (unrooted, non-trivial splits):
 * @code
 *   int rf = t1.robinson_foulds(t2);
 * @endcode
 */
class PhyloTree {
public:
	std::shared_ptr<TreeNode> root;

	PhyloTree() : root(std::make_shared<TreeNode>()) {}
	explicit PhyloTree(std::shared_ptr<TreeNode> r) : root(std::move(r)) {}

	/** Parse a Newick string and return the corresponding PhyloTree.
	 *  Throws std::invalid_argument if the string contains trailing content
	 *  after the optional terminating semicolon.
	 */
	static PhyloTree from_newick(const std::string& newick) {
		size_t pos = 0;
		auto r = detail::parse_node(newick, pos);
		detail::skip_ws(newick, pos);
		if (pos < newick.size() && newick[pos] == ';') ++pos;
		detail::skip_ws(newick, pos);
		if (pos < newick.size())
			throw std::invalid_argument(
				"Newick parse error: unexpected trailing content: '"
				+ newick.substr(pos) + "'");
		return PhyloTree(std::move(r));
	}

	/** Serialise the tree as a Newick string (trailing semicolon included). */
	std::string to_newick() const {
		return root ? root->to_newick() : ";";
	}

	/** Names of all leaves in postorder. */
	std::vector<std::string> get_leaf_names() const {
		return root ? root->get_leaf_names() : std::vector<std::string>{};
	}

	/** Prune the tree to keep only the listed taxa (see TreeNode::prune). */
	void prune(const std::vector<std::string>& taxa_to_keep) {
		if (root) root->prune(taxa_to_keep);
	}

	/**
	 * Return all rooted triplets displayed by this tree.
	 *
	 * Each triplet is represented as {a, b, c} where a and b form the cherry,
	 * i.e. the triplet topology is ((a,b),c).  The order of a and b within
	 * the cherry is unspecified.
	 *
	 * Algorithm: for each internal node v, and for each child ci of v whose
	 * subtree contains at least two leaves, every pair (a,b) within ci's
	 * subtree forms a cherry with every leaf c in any other child cj of v.
	 * This is correct because lca(a,b,c) = v and lca(a,b) is strictly below
	 * v, and each such triplet is generated exactly once.
	 */
	std::vector<std::array<std::string, 3>> get_triplets() const {
		std::vector<std::array<std::string, 3>> triplets;
		if (!root) return triplets;

		for (auto* node : root->traverse()) {
			if (node->is_leaf()) continue;

			// Collect leaf names for each child's subtree.
			std::vector<std::vector<std::string>> child_leaves;
			child_leaves.reserve(node->children.size());
			for (const auto& child : node->children)
				child_leaves.push_back(child->get_leaf_names());

			for (size_t i = 0; i < child_leaves.size(); ++i) {
				const auto& cherry = child_leaves[i];
				if (cherry.size() < 2) continue;

				// Build the flat list of outgroup leaves from all other children.
				std::vector<const std::string*> outgroup;
				for (size_t j = 0; j < child_leaves.size(); ++j) {
					if (j == i) continue;
					for (const auto& name : child_leaves[j])
						outgroup.push_back(&name);
				}

				for (size_t a = 0; a < cherry.size(); ++a)
					for (size_t b = a + 1; b < cherry.size(); ++b)
						for (const auto* c : outgroup)
							triplets.push_back({cherry[a], cherry[b], *c});
			}
		}
		return triplets;
	}

	/**
	 * Unrooted Robinson-Foulds distance between this tree and `other`.
	 *
	 * Returns the number of non-trivial bipartitions (splits) that differ
	 * between the two trees — i.e. |splits_only_here| + |splits_only_there|.
	 * Trivial splits (isolating a single leaf) are excluded.
	 */
	int robinson_foulds(const PhyloTree& other) const {
		auto s1 = bipartitions();
		auto s2 = other.bipartitions();
		int d = 0;
		for (const auto& s : s1) if (!s2.count(s)) ++d;
		for (const auto& s : s2) if (!s1.count(s)) ++d;
		return d;
	}

private:
	/**
	 * Compute the set of non-trivial bipartitions of this tree.
	 *
	 * Each bipartition is represented by the lexicographically smaller of its
	 * two halves, ensuring each split has exactly one canonical key regardless
	 * of which internal node we encounter it from.
	 */
	std::set<std::set<std::string>> bipartitions() const {
		if (!root) return {};

		auto all_leaves = root->get_leaf_names();
		std::set<std::string> all_set(all_leaves.begin(), all_leaves.end());
		std::set<std::set<std::string>> splits;

		for (auto* node : root->traverse()) {
			if (node->is_leaf()) continue;
			if (node == root.get()) continue;

			auto sub = node->get_leaf_names();
			std::set<std::string> side1(sub.begin(), sub.end());

			// Skip trivial bipartitions.
			if (side1.size() <= 1 || side1.size() >= all_set.size() - 1) continue;

			std::set<std::string> side2;
			for (const auto& l : all_set)
				if (!side1.count(l)) side2.insert(l);

			// Canonicalize by lexicographic order so both internal nodes on
			// either side of the same split produce the same key.
			splits.insert(side1 < side2 ? side1 : side2);
		}
		return splits;
	}
};

} // namespace phylotriplet
