/**
 * newick2triplets — enumerate rooted triplets from a tree and print them in
 * the weighted-trees file format accepted by the phylotriplet inference
 * algorithm ("newick_triplet weight" per line).
 *
 * Usage (three input modes):
 *
 *   newick2triplets "((a,b),c);"          # single Newick via command-line arg
 *   echo "((a,b),c);" | newick2triplets   # single Newick via stdin
 *   newick2triplets --file trees.txt      # weighted-trees file (Newick + weight)
 *
 * In the single-Newick modes the output weight for every triplet is 1.0.
 * In file mode each triplet inherits the weight of its source tree.
 *
 * Output format (one triplet per line):
 *   ((a,b),c); <weight>
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include "phylotriplet/file_io.hpp"
#include "phylotriplet/triplet_estimator.hpp"

using namespace phylotriplet;

static void print_triplets(const std::vector<WeightedTriplet>& wts)
{
	for (const auto& wt : wts) {
		const auto& t = wt.triplet;
		// t encodes a,b|c  →  topology ((a,b),c)
		std::cout << "((" << t[0] << "," << t[1] << ")," << t[2] << "); "
		          << wt.weight << "\n";
	}
}

int main(int argc, char* argv[])
{
	try {
		if (argc == 3 && std::string(argv[1]) == "--file") {
			// Weighted-trees file mode
			std::string path(argv[2]);
			auto trees = read_weighted_trees_from_file(path);
			auto wts   = extract_weighted_triplets(trees);
			print_triplets(wts);

		} else if (argc == 2 && std::string(argv[1]) != "--file") {
			// Single Newick from command-line argument
			std::string newick(argv[1]);
			PhyloTree tree = PhyloTree::from_newick(newick);
			std::vector<WeightedTriplet> wts;
			for (const auto& t : tree.get_triplets())
				wts.emplace_back(t, 1.0);
			print_triplets(wts);

		} else if (argc == 1) {
			// Single Newick from stdin
			std::string newick;
			std::getline(std::cin, newick);
			if (newick.empty()) {
				std::cerr << "Error: no input on stdin\n";
				return 1;
			}
			PhyloTree tree = PhyloTree::from_newick(newick);
			std::vector<WeightedTriplet> wts;
			for (const auto& t : tree.get_triplets())
				wts.emplace_back(t, 1.0);
			print_triplets(wts);

		} else {
			std::cerr <<
				"Usage:\n"
				"  newick2triplets \"<newick>\"       single Newick from argument\n"
				"  newick2triplets --file <path>    weighted-trees file\n"
				"  echo \"<newick>\" | newick2triplets  single Newick from stdin\n";
			return 1;
		}

	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}
	return 0;
}