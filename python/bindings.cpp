#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "phylotriplet/triplet_estimator.hpp"
#include "phylotriplet/file_io.hpp"

namespace py = pybind11;
using namespace phylotriplet;

PYBIND11_MODULE(phylotriplet, m) {
    m.doc() = "Phylogenetic tree estimation from weighted input trees";

    py::class_<PhyloTree>(m, "PhyloTree")
        .def_static("from_newick", &PhyloTree::from_newick, py::arg("newick"),
                    "Parse a Newick string and return a PhyloTree.")
        .def("to_newick", &PhyloTree::to_newick,
             "Serialise the tree to a Newick string.")
        .def("get_leaf_names", &PhyloTree::get_leaf_names,
             "Return leaf names in postorder.")
        .def("get_triplets", &PhyloTree::get_triplets,
             "Return all rooted triplets as [a, b, c] lists where a,b|c "
             "(i.e. the topology is ((a,b),c)).")
        .def("robinson_foulds", &PhyloTree::robinson_foulds, py::arg("other"),
             "Unrooted Robinson-Foulds distance to another PhyloTree.")
        .def("__repr__", [](const PhyloTree& t) {
            return "PhyloTree('" + t.to_newick() + "')";
        });

    py::class_<WeightedTree>(m, "WeightedTree")
        .def(py::init<const std::string&, double>(),
             py::arg("newick"), py::arg("weight") = 1.0,
             "Construct from a Newick string and optional weight.")
        .def(py::init<PhyloTree, double>(),
             py::arg("tree"), py::arg("weight") = 1.0,
             "Construct from a PhyloTree and optional weight.")
        .def_property("newick",
            [](const WeightedTree& wt) { return wt.tree.to_newick(); },
            [](WeightedTree& wt, const std::string& nwk) {
                wt.tree = PhyloTree::from_newick(nwk);
            },
            "Newick string representation of the stored tree.")
        .def_property("tree",
            [](const WeightedTree& wt) { return wt.tree; },
            [](WeightedTree& wt, const PhyloTree& t) { wt.tree = t; },
            "The stored PhyloTree object.")
        .def_readwrite("weight", &WeightedTree::weight)
        .def("__repr__", [](const WeightedTree& wt) {
            return "WeightedTree(newick=\"" + wt.tree.to_newick()
                   + "\", weight=" + std::to_string(wt.weight) + ")";
        });

    // I/O utility functions
    m.def("read_weighted_trees_from_file", &read_weighted_trees_from_file,
          py::arg("filepath"),
          "Read weighted trees from a file.\n\n"
          "Format: one 'newick_tree weight' per line.\n\n"
          "Args:\n"
          "    filepath: Path to input file\n\n"
          "Returns:\n"
          "    list of WeightedTree");

    m.def("read_weighted_trees_from_string", &read_weighted_trees_from_string,
          py::arg("content"),
          "Parse weighted trees from a string.\n\n"
          "Format: one 'newick_tree weight' per line.\n\n"
          "Args:\n"
          "    content: String with weighted trees\n\n"
          "Returns:\n"
          "    list of WeightedTree");

    // Tree estimation
    m.def("estimate_tree_from_file", &estimate_tree_from_file,
          py::arg("filepath"),
          "Estimate a phylogenetic tree from a weighted trees file.\n\n"
          "Args:\n"
          "    filepath: Path to input file (one 'newick_tree weight' per line)\n\n"
          "Returns:\n"
          "    str: Newick string of the estimated tree");

    m.def("estimate_tree_from_string", &estimate_tree_from_string,
          py::arg("content"),
          "Estimate a phylogenetic tree from a weighted trees string.\n\n"
          "Args:\n"
          "    content: String with weighted trees (one 'newick_tree weight' per line)\n\n"
          "Returns:\n"
          "    str: Newick string of the estimated tree");

}