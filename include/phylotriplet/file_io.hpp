#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tree.hpp"
#include "warn.hpp"

namespace phylotriplet {

/**
 * A phylogenetic tree with an associated weight.
 * Initialise from a Newick string; the string is parsed once on construction.
 */
struct WeightedTree {
  PhyloTree tree;  // parsed phylogenetic tree
  double weight;   // tree weight

  WeightedTree(const std::string& nwk, double w = 1.0)
      : tree(PhyloTree::from_newick(nwk)), weight(w) {}

  explicit WeightedTree(PhyloTree t, double w = 1.0)
      : tree(std::move(t)), weight(w) {}

  /** Serialise the stored tree back to a Newick string. */
  std::string newick() const { return tree.to_newick(); }
};

namespace detail {

/**
 * Parse weighted trees from a stream.
 * @param stream  Input stream to read from.
 * @param source  Human-readable name used in error messages (e.g. file path
 *                or "<string>").
 */
inline std::vector<WeightedTree> parse_weighted_trees_from_stream(
    std::istream& stream, const std::string& source) {
  std::vector<WeightedTree> trees;
  std::string line;
  int line_num = 0;
  bool weight_warning_already_displayed = false;

  while (std::getline(stream, line)) {
    line_num++;

    // Normalize whitespace: collapse multiple spaces/tabs into single space
    // and strip leading whitespace.
    std::string clean_line;
    clean_line.reserve(line.size());
    bool prev_was_whitespace = true;  // strip leading whitespace
    for (char c : line) {
      if (std::isspace(static_cast<unsigned char>(c))) {
        if (!prev_was_whitespace) {
          clean_line += ' ';
          prev_was_whitespace = true;
        }
      } else {
        clean_line += c;
        prev_was_whitespace = false;
      }
    }
    line = std::move(clean_line);

    // Strip any trailing space
    if (!line.empty() &&
        std::isspace(static_cast<unsigned char>(line.back()))) {
      line.pop_back();
    }

    // Skip empty lines and comments.
    if (line.empty() || line[0] == '#') {
      continue;
    }

    double weight = 1.0;
    // Split on last space for weight.
    size_t space_separator = line.find_last_of(' ');
    if (space_separator != std::string::npos) {
      std::string weight_str = line.substr(space_separator);

      // Parse weight.
      try {
        weight = std::stod(weight_str);
      } catch (...) {
        throw std::runtime_error("Invalid weight at line " +
                                 std::to_string(line_num) + " of " + source +
                                 ": '" + weight_str + "'");
      }
    } else {
      if (!weight_warning_already_displayed) {
        PHYLOTRIPLET_WARN("Warning: no weight found at line " +
                          std::to_string(line_num) + " of " + source +
                          ", defaulting to 1.0.");
      }
    }

    std::string newick = line.substr(0, space_separator);

    // Trim any remaining whitespace from the newick string.
    newick.erase(std::remove(newick.begin(), newick.end(), ' '), newick.end());

    // Basic Newick validation: must end with semicolon.
    if (newick.empty()) {
      throw std::runtime_error("Empty Newick string at line " +
                               std::to_string(line_num) + " of " + source);
    }
    if (newick.back() != ';') {
      newick += ';';
      PHYLOTRIPLET_WARN("Warning: Newick string missing semicolon at line " +
                        std::to_string(line_num) + " of " + source + ".");
    }

    try {
      trees.emplace_back(newick, weight);
    } catch (const std::exception& e) {
      throw std::runtime_error("Error on line " + std::to_string(line_num) +
                               " of " + source + ": " + e.what());
    }
  }

  return trees;
}

}  // namespace detail

/**
 * Read weighted trees from file.
 * Format: newick_tree weight (one per line)
 *
 * @param filepath Path to the input file
 * @return Vector of WeightedTree objects
 * @throws std::runtime_error if file cannot be opened, format is invalid,
 *         or the file contains no trees
 */
inline std::vector<WeightedTree> read_weighted_trees_from_file(
    const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + filepath);
  }

  auto trees = detail::parse_weighted_trees_from_stream(file, filepath);

  if (trees.empty()) {
    throw std::runtime_error("No trees found in file: " + filepath);
  }

  return trees;
}

/**
 * Parse weighted trees from string content.
 * Format: newick_tree weight (one per line)
 *
 * @param content String containing weighted trees
 * @return Vector of WeightedTree objects
 * @throws std::runtime_error if format is invalid
 */
inline std::vector<WeightedTree> read_weighted_trees_from_string(
    const std::string& content) {
  std::istringstream stream(content);
  return detail::parse_weighted_trees_from_stream(stream, "<string>");
}

}  // namespace phylotriplet