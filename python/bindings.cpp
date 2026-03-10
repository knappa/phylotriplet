#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define PHYLOTRIPLET_PY_WARNINGS
#include "phylotriplet/file_io.hpp"
#include "phylotriplet/sequence_analysis.hpp"
#include "phylotriplet/sequence_io.hpp"
#include "phylotriplet/triplet_estimator.hpp"

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
      .def(py::init<const std::string&, double>(), py::arg("newick"),
           py::arg("weight") = 1.0,
           "Construct from a Newick string and optional weight.")
      .def(py::init<PhyloTree, double>(), py::arg("tree"),
           py::arg("weight") = 1.0,
           "Construct from a PhyloTree and optional weight.")
      .def_property(
          "newick", [](const WeightedTree& wt) { return wt.tree.to_newick(); },
          [](WeightedTree& wt, const std::string& nwk) {
            wt.tree = PhyloTree::from_newick(nwk);
          },
          "Newick string representation of the stored tree.")
      .def_property(
          "tree", [](const WeightedTree& wt) { return wt.tree; },
          [](WeightedTree& wt, const PhyloTree& t) { wt.tree = t; },
          "The stored PhyloTree object.")
      .def_readwrite("weight", &WeightedTree::weight)
      .def("__repr__", [](const WeightedTree& wt) {
        return "WeightedTree(newick=\"" + wt.tree.to_newick() +
               "\", weight=" + std::to_string(wt.weight) + ")";
      });

  // -- Sequence alignment ----------------------------------------------------

  py::class_<TuplePatterns>(
      m, "TuplePatterns",
      "Site-pattern counts for an ordered k-tuple of taxa (k >= 2).\n\n"
      "taxa         -- sorted list of k taxon names\n"
      "site_pattern -- list of 4^k floats.\n"
      "               site_pattern[b0*4^(k-1) + b1*4^(k-2) + ... + b_{k-1}]\n"
      "               is the count for taxa[0]=b0, ..., taxa[k-1]=b_{k-1},\n"
      "               with 0=A 1=C 2=G 3=T\n"
      "as_array()   -- site_pattern as a numpy array of shape (4,)*k\n"
      "n_total      -- total alignment columns (including gap columns)\n"
      "n_counted    -- columns with no gap in any taxon")
      .def_readonly("taxa", &TuplePatterns::taxa, "Sorted list of taxon names.")
      .def_readonly("site_pattern", &TuplePatterns::site_pattern,
                    "4^k site-pattern counts.")
      .def_property_readonly("order", &TuplePatterns::order,
                             "Number of taxa in the tuple (k).")
      .def_readonly("n_total", &TuplePatterns::n_total,
                    "Total alignment columns.")
      .def_readonly("n_counted", &TuplePatterns::n_counted,
                    "Columns with no gap in any taxon.")
      .def(
          "as_array",
          [](const TuplePatterns& p) {
            const size_t k = p.order();
            std::vector<py::ssize_t> shape(k, 4);
            py::array_t<double> arr(shape);
            auto raw = arr.mutable_data();
            for (size_t i = 0; i < p.site_pattern.size(); ++i)
              raw[i] = p.site_pattern[i];
            return arr;
          },
          "Return site-pattern counts as a numpy array of shape (4,)*k.\n"
          "Requires numpy (install with: pip install numpy).")
      .def("__repr__", [](const TuplePatterns& p) {
        std::string s = "TuplePatterns(taxa=[";
        for (size_t i = 0; i < p.taxa.size(); ++i) {
          if (i > 0) s += ", ";
          s += "'" + p.taxa[i] + "'";
        }
        s += "], n_counted=" + std::to_string(p.n_counted) + ")";
        return s;
      });
  // Backward-compatible alias
  m.attr("TripletSitePatterns") = m.attr("TuplePatterns");

  py::class_<GeneAlignment>(m, "GeneAlignment",
                            "A single-gene alignment of named sequences.\n\n"
                            "gene_name    -- name of this gene\n"
                            "num_taxa     -- number of sequences\n"
                            "length       -- alignment length in columns\n"
                            "taxon_names() -- sorted list of taxon names")
      .def_readonly("gene_name", &GeneAlignment::gene_name,
                    "Name of this gene.")
      .def_property_readonly("num_taxa", &GeneAlignment::num_taxa,
                             "Number of sequences.")
      .def_property_readonly("length", &GeneAlignment::length,
                             "Alignment length in columns.")
      .def("taxon_names", &GeneAlignment::taxon_names,
           "List of taxon names (in sequence order).")
      .def(
          "triplet_site_patterns",
          [](const GeneAlignment& gene) { return triplet_site_patterns(gene); },
          "All C(n,3) triplet site-pattern counts for this gene.")
      .def(
          "triplet_site_pattern",
          [](const GeneAlignment& gene, const std::string& t1,
             const std::string& t2, const std::string& t3) {
            return triplet_site_pattern(gene, t1, t2, t3);
          },
          py::arg("taxon1"), py::arg("taxon2"), py::arg("taxon3"),
          "Site-pattern counts for a specific triplet.\n\n"
          "Raises ValueError if any taxon is absent.")
      .def(
          "quartet_site_patterns",
          [](const GeneAlignment& gene) { return quartet_site_patterns(gene); },
          "All C(n,4) quartet site-pattern counts for this gene.")
      .def(
          "quartet_site_pattern",
          [](const GeneAlignment& gene, const std::string& t1,
             const std::string& t2, const std::string& t3,
             const std::string& t4) {
            return quartet_site_pattern(gene, t1, t2, t3, t4);
          },
          py::arg("taxon1"), py::arg("taxon2"), py::arg("taxon3"),
          py::arg("taxon4"),
          "Site-pattern counts for a specific quartet.\n\n"
          "Raises ValueError if any taxon is absent.")
      .def(
          "site_pattern",
          [](const GeneAlignment& gene, const std::vector<std::string>& taxa) {
            return site_pattern(gene, taxa);
          },
          py::arg("taxa"),
          "Site-pattern counts for a specific k-tuple of taxa (k >= 2).\n\n"
          "Raises ValueError if any taxon is absent.")
      .def("__repr__", [](const GeneAlignment& g) {
        return "GeneAlignment(gene_name='" + g.gene_name +
               "', num_taxa=" + std::to_string(g.num_taxa()) +
               ", length=" + std::to_string(g.length()) + ")";
      });

  py::class_<MultiGeneAlignment>(
      m, "MultiGeneAlignment",
      "A collection of per-gene alignments.\n\n"
      "Taxa may differ between genes (partial data is allowed).")
      .def_property_readonly("num_genes", &MultiGeneAlignment::num_genes,
                             "Number of gene alignments.")
      .def("gene_names", &MultiGeneAlignment::gene_names,
           "List of gene names in order.")
      .def("all_taxon_names", &MultiGeneAlignment::all_taxon_names,
           "Sorted list of all taxon names across all genes.")
      .def("__len__", &MultiGeneAlignment::num_genes)
      .def(
          "__getitem__",
          [](const MultiGeneAlignment& aln,
             py::object key) -> const GeneAlignment& {
            if (py::isinstance<py::int_>(key)) {
              long i = key.cast<long>();
              if (i < 0) i += static_cast<long>(aln.num_genes());
              if (i < 0 || static_cast<size_t>(i) >= aln.num_genes())
                throw py::index_error("gene index out of range");
              return aln.genes[static_cast<size_t>(i)];
            }
            if (py::isinstance<py::str>(key)) {
              const std::string name = key.cast<std::string>();
              for (const auto& g : aln.genes)
                if (g.gene_name == name) return g;
              throw py::key_error("gene '" + name + "' not found");
            }
            throw py::type_error("gene index must be int or str");
          },
          py::return_value_policy::reference_internal,
          "Access a gene by integer index or name string.")
      .def(
          "__iter__",
          [](const MultiGeneAlignment& aln) {
            return py::make_iterator(aln.genes.begin(), aln.genes.end());
          },
          py::keep_alive<0, 1>(), "Iterate over GeneAlignment objects.")
      .def(
          "triplet_site_patterns",
          [](const MultiGeneAlignment& aln) {
            return triplet_site_patterns(aln);
          },
          "Compute aggregate site-pattern counts for all C(n,3) taxon "
          "triplets across all genes.\n\n"
          "Returns a list of TuplePatterns, one per triplet.  Contributions\n"
          "from a gene are included only if it contains all three taxa;\n"
          "n_total and n_counted are summed accordingly.")
      .def(
          "triplet_site_pattern",
          [](const MultiGeneAlignment& aln, const std::string& t1,
             const std::string& t2, const std::string& t3) {
            return triplet_site_pattern(aln, t1, t2, t3);
          },
          py::arg("taxon1"), py::arg("taxon2"), py::arg("taxon3"),
          "Compute aggregate site-pattern counts for a specific taxon "
          "triplet.\n\n"
          "Taxa names may be given in any order.  Genes that do not contain\n"
          "all three taxa are skipped with a warning.\n\n"
          "Args:\n"
          "    taxon1, taxon2, taxon3: taxon names (any order)\n\n"
          "Returns:\n"
          "    TuplePatterns")
      .def(
          "quartet_site_patterns",
          [](const MultiGeneAlignment& aln) {
            return quartet_site_patterns(aln);
          },
          "Compute aggregate site-pattern counts for all C(n,4) taxon "
          "quartets across all genes.\n\n"
          "Returns a list of TuplePatterns, one per quartet.  Only genes\n"
          "that contain all four taxa contribute.")
      .def(
          "quartet_site_pattern",
          [](const MultiGeneAlignment& aln, const std::string& t1,
             const std::string& t2, const std::string& t3,
             const std::string& t4) {
            return quartet_site_pattern(aln, t1, t2, t3, t4);
          },
          py::arg("taxon1"), py::arg("taxon2"), py::arg("taxon3"),
          py::arg("taxon4"),
          "Compute aggregate site-pattern counts for a specific taxon "
          "quartet.\n\n"
          "Taxa names may be given in any order.  Genes missing any taxon\n"
          "are skipped with a warning.\n\n"
          "Args:\n"
          "    taxon1, taxon2, taxon3, taxon4: taxon names (any order)\n\n"
          "Returns:\n"
          "    TuplePatterns")
      .def(
          "site_pattern",
          [](const MultiGeneAlignment& aln,
             const std::vector<std::string>& taxa) {
            return site_pattern(aln, taxa);
          },
          py::arg("taxa"),
          "Compute aggregate site-pattern counts for a k-tuple of taxa "
          "(k >= 2).\n\n"
          "Taxa names may be given in any order.  Genes missing any taxon\n"
          "are skipped with a warning.\n\n"
          "Args:\n"
          "    taxa: list of taxon names (any order, length >= 2)\n\n"
          "Returns:\n"
          "    TuplePatterns")
      .def("__repr__", [](const MultiGeneAlignment& aln) {
        return "MultiGeneAlignment(num_genes=" +
               std::to_string(aln.num_genes()) +
               ", num_taxa=" + std::to_string(aln.all_taxon_names().size()) +
               ")";
      });

  // -- I/O: sequence alignments ----------------------------------------------

  m.def("read_nexus_file", &read_nexus_file, py::arg("filepath"),
        "Read a NEXUS alignment file.\n\n"
        "Supports multiple DATA/CHARACTERS blocks (one per gene), a single\n"
        "DATA block split by a SETS block with CHARSET commands, sequential\n"
        "and interleaved MATRIX formats, and DNA/RNA datatypes.\n\n"
        "Args:\n"
        "    filepath: Path to the NEXUS file\n\n"
        "Returns:\n"
        "    MultiGeneAlignment");

  m.def(
      "read_phylip_file",
      [](const std::string& filepath, const std::string& gene_name,
         bool interleaved) {
        MultiGeneAlignment aln;
        aln.genes.push_back(read_phylip_file(filepath, gene_name, interleaved));
        return aln;
      },
      py::arg("filepath"), py::arg("gene_name") = std::string(""),
      py::arg("interleaved") = false,
      "Read a single PHYLIP alignment file.\n\n"
      "Args:\n"
      "    filepath:   Path to the PHYLIP file\n"
      "    gene_name:  Name for this gene; defaults to the filename stem\n"
      "    interleaved: Set True for interleaved PHYLIP format\n\n"
      "Returns:\n"
      "    MultiGeneAlignment (with one gene)");

  m.def("read_phylip_files", &read_phylip_files, py::arg("filepaths"),
        py::arg("gene_names") = std::vector<std::string>{},
        py::arg("interleaved") = false,
        "Read multiple PHYLIP alignment files as a multi-gene alignment.\n\n"
        "Args:\n"
        "    filepaths:   List of paths, one per gene\n"
        "    gene_names:  Optional list of gene names; defaults to filename "
        "stems\n"
        "    interleaved: Applied uniformly to all files\n\n"
        "Returns:\n"
        "    MultiGeneAlignment");

  // -- I/O: weighted trees ---------------------------------------------------

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

  // -- Tree estimation -------------------------------------------------------

  m.def(
      "estimate_tree_from_file", &estimate_tree_from_file, py::arg("filepath"),
      "Estimate a phylogenetic tree from a weighted trees file.\n\n"
      "Args:\n"
      "    filepath: Path to input file (one 'newick_tree weight' per line)\n\n"
      "Returns:\n"
      "    str: Newick string of the estimated tree");

  m.def("estimate_tree_from_string", &estimate_tree_from_string,
        py::arg("content"),
        "Estimate a phylogenetic tree from a weighted trees string.\n\n"
        "Args:\n"
        "    content: String with weighted trees (one 'newick_tree weight' per "
        "line)\n\n"
        "Returns:\n"
        "    str: Newick string of the estimated tree");
}