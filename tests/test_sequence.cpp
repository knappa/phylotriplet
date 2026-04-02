#define CATCH_CONFIG_MAIN
#include <sstream>

#include "catch.hpp"
#include "phylotriplet/sequence.hpp"
#include "phylotriplet/sequence_analysis.hpp"
#include "phylotriplet/sequence_io.hpp"

using namespace phylotriplet;

// --- helpers -----------------------------------------------------------------

// Build a GeneAlignment directly from (taxon_name, bases_string) pairs.
static GeneAlignment make_gene(
    const std::vector<std::pair<std::string, std::string>>& taxa,
    const std::string& gene_name = "gene") {
  GeneAlignment gene;
  gene.gene_name = gene_name;
  for (const auto& kv : taxa) {
    Sequence s;
    s.taxon = kv.first;
    for (char c : kv.second) s.bases.push_back(char_to_nucleotide(c));
    gene.sequences.push_back(std::move(s));
  }
  return gene;
}

// Sum all 64 site-pattern counts.
static double total_pattern_mass(const TuplePatterns& pats) {
  double t = 0.0;
  for (int k = 0; k < 64; ++k) t += pats.site_pattern[k];
  return t;
}

// --- unambiguous_base_index --------------------------------------------------

TEST_CASE("unambiguous_base_index: A=0, C=1, G=2, T=3") {
  REQUIRE(unambiguous_base_index(NUC_A) == 0);
  REQUIRE(unambiguous_base_index(NUC_C) == 1);
  REQUIRE(unambiguous_base_index(NUC_G) == 2);
  REQUIRE(unambiguous_base_index(NUC_T) == 3);
}

TEST_CASE("unambiguous_base_index: gap returns -1") {
  REQUIRE(unambiguous_base_index(NUC_GAP) == -1);
}

TEST_CASE("unambiguous_base_index: ambiguous codes return -1") {
  REQUIRE(unambiguous_base_index(NUC_R) == -1);
  REQUIRE(unambiguous_base_index(NUC_Y) == -1);
  REQUIRE(unambiguous_base_index(NUC_S) == -1);
  REQUIRE(unambiguous_base_index(NUC_W) == -1);
  REQUIRE(unambiguous_base_index(NUC_K) == -1);
  REQUIRE(unambiguous_base_index(NUC_M) == -1);
  REQUIRE(unambiguous_base_index(NUC_B) == -1);
  REQUIRE(unambiguous_base_index(NUC_D) == -1);
  REQUIRE(unambiguous_base_index(NUC_H) == -1);
  REQUIRE(unambiguous_base_index(NUC_V) == -1);
  REQUIRE(unambiguous_base_index(NUC_N) == -1);
}

// --- char_to_nucleotide ------------------------------------------------------

TEST_CASE("char_to_nucleotide: unambiguous bases (upper and lower case)") {
  REQUIRE(char_to_nucleotide('A') == NUC_A);
  REQUIRE(char_to_nucleotide('C') == NUC_C);
  REQUIRE(char_to_nucleotide('G') == NUC_G);
  REQUIRE(char_to_nucleotide('T') == NUC_T);
  REQUIRE(char_to_nucleotide('a') == NUC_A);
  REQUIRE(char_to_nucleotide('c') == NUC_C);
  REQUIRE(char_to_nucleotide('g') == NUC_G);
  REQUIRE(char_to_nucleotide('t') == NUC_T);
}

TEST_CASE("char_to_nucleotide: U and u map to T") {
  REQUIRE(char_to_nucleotide('U') == NUC_T);
  REQUIRE(char_to_nucleotide('u') == NUC_T);
}

TEST_CASE("char_to_nucleotide: IUPAC ambiguity codes") {
  REQUIRE(char_to_nucleotide('R') == (NUC_A | NUC_G));
  REQUIRE(char_to_nucleotide('Y') == (NUC_C | NUC_T));
  REQUIRE(char_to_nucleotide('S') == (NUC_C | NUC_G));
  REQUIRE(char_to_nucleotide('W') == (NUC_A | NUC_T));
  REQUIRE(char_to_nucleotide('K') == (NUC_G | NUC_T));
  REQUIRE(char_to_nucleotide('M') == (NUC_A | NUC_C));
  REQUIRE(char_to_nucleotide('B') == (NUC_C | NUC_G | NUC_T));
  REQUIRE(char_to_nucleotide('D') == (NUC_A | NUC_G | NUC_T));
  REQUIRE(char_to_nucleotide('H') == (NUC_A | NUC_C | NUC_T));
  REQUIRE(char_to_nucleotide('V') == (NUC_A | NUC_C | NUC_G));
  REQUIRE(char_to_nucleotide('N') == (NUC_A | NUC_C | NUC_G | NUC_T));
  REQUIRE(char_to_nucleotide('X') == (NUC_A | NUC_C | NUC_G | NUC_T));
}

TEST_CASE("char_to_nucleotide: IUPAC codes are case-insensitive") {
  REQUIRE(char_to_nucleotide('r') == char_to_nucleotide('R'));
  REQUIRE(char_to_nucleotide('n') == char_to_nucleotide('N'));
}

TEST_CASE("char_to_nucleotide: gap and unknown characters return NUC_GAP") {
  REQUIRE(char_to_nucleotide('-') == NUC_GAP);
  REQUIRE(char_to_nucleotide('?') == NUC_GAP);
  REQUIRE(char_to_nucleotide('Z') == NUC_GAP);
  REQUIRE(char_to_nucleotide(' ') == NUC_GAP);
}

// --- nucleotide_distribution -------------------------------------------------

TEST_CASE("nucleotide_distribution: unambiguous bases have probability 1") {
  NucleotideDistribution d;

  d = nucleotide_distribution(NUC_A);
  REQUIRE(d[0] == Approx(1.0));
  REQUIRE(d[1] == Approx(0.0));
  REQUIRE(d[2] == Approx(0.0));
  REQUIRE(d[3] == Approx(0.0));

  d = nucleotide_distribution(NUC_C);
  REQUIRE(d[0] == Approx(0.0));
  REQUIRE(d[1] == Approx(1.0));
  REQUIRE(d[2] == Approx(0.0));
  REQUIRE(d[3] == Approx(0.0));

  d = nucleotide_distribution(NUC_G);
  REQUIRE(d[0] == Approx(0.0));
  REQUIRE(d[1] == Approx(0.0));
  REQUIRE(d[2] == Approx(1.0));
  REQUIRE(d[3] == Approx(0.0));

  d = nucleotide_distribution(NUC_T);
  REQUIRE(d[0] == Approx(0.0));
  REQUIRE(d[1] == Approx(0.0));
  REQUIRE(d[2] == Approx(0.0));
  REQUIRE(d[3] == Approx(1.0));
}

TEST_CASE("nucleotide_distribution: N is uniform over all four bases") {
  NucleotideDistribution d = nucleotide_distribution(NUC_N);
  REQUIRE(d[0] == Approx(0.25));
  REQUIRE(d[1] == Approx(0.25));
  REQUIRE(d[2] == Approx(0.25));
  REQUIRE(d[3] == Approx(0.25));
}

TEST_CASE("nucleotide_distribution: R splits equally between A and G") {
  NucleotideDistribution d = nucleotide_distribution(NUC_R);  // A|G
  REQUIRE(d[0] == Approx(0.5));                               // A
  REQUIRE(d[1] == Approx(0.0));                               // C
  REQUIRE(d[2] == Approx(0.5));                               // G
  REQUIRE(d[3] == Approx(0.0));                               // T
}

TEST_CASE("nucleotide_distribution: GAP returns all zeros") {
  NucleotideDistribution d = nucleotide_distribution(NUC_GAP);
  REQUIRE(d[0] == Approx(0.0));
  REQUIRE(d[1] == Approx(0.0));
  REQUIRE(d[2] == Approx(0.0));
  REQUIRE(d[3] == Approx(0.0));
}

TEST_CASE("nucleotide_distribution: non-gap distribution always sums to 1") {
  const uint8_t codes[] = {NUC_A, NUC_C, NUC_G, NUC_T, NUC_R,
                           NUC_Y, NUC_S, NUC_W, NUC_K, NUC_M,
                           NUC_B, NUC_D, NUC_H, NUC_V, NUC_N};
  for (uint8_t b : codes) {
    NucleotideDistribution d = nucleotide_distribution(b);
    double s = d[0] + d[1] + d[2] + d[3];
    REQUIRE(s == Approx(1.0));
  }
}

// --- read_phylip -------------------------------------------------------------

TEST_CASE("read_phylip sequential: taxon names, bases, gene_name") {
  const std::string phy =
      "3 4\n"
      "Alpha ACGT\n"
      "Beta  TTAG\n"
      "Gamma CGCA\n";
  std::istringstream ss(phy);
  GeneAlignment gene = read_phylip(ss, "mygene");

  REQUIRE(gene.gene_name == "mygene");
  REQUIRE(gene.num_taxa() == 3);
  REQUIRE(gene.length() == 4);
  REQUIRE(gene.sequences[0].taxon == "Alpha");
  REQUIRE(gene.sequences[1].taxon == "Beta");
  REQUIRE(gene.sequences[2].taxon == "Gamma");

  // Alpha: A C G T
  REQUIRE(gene.sequences[0].bases[0] == NUC_A);
  REQUIRE(gene.sequences[0].bases[1] == NUC_C);
  REQUIRE(gene.sequences[0].bases[2] == NUC_G);
  REQUIRE(gene.sequences[0].bases[3] == NUC_T);

  // Beta: T T A G
  REQUIRE(gene.sequences[1].bases[0] == NUC_T);
  REQUIRE(gene.sequences[1].bases[2] == NUC_A);
}

TEST_CASE("read_phylip sequential: IUPAC ambiguity codes preserved") {
  const std::string phy =
      "2 4\n"
      "T1 ACNR\n"
      "T2 TYWK\n";
  std::istringstream ss(phy);
  GeneAlignment gene = read_phylip(ss);

  REQUIRE(gene.sequences[0].bases[2] == NUC_N);
  REQUIRE(gene.sequences[0].bases[3] == NUC_R);
  REQUIRE(gene.sequences[1].bases[1] == NUC_Y);
  REQUIRE(gene.sequences[1].bases[3] == NUC_K);
}

TEST_CASE("read_phylip sequential: gap characters encoded as NUC_GAP") {
  const std::string phy =
      "2 3\n"
      "T1 A-G\n"
      "T2 C?T\n";
  std::istringstream ss(phy);
  GeneAlignment gene = read_phylip(ss);

  REQUIRE(gene.sequences[0].bases[1] == NUC_GAP);
  REQUIRE(gene.sequences[1].bases[1] == NUC_GAP);
}

TEST_CASE("read_phylip interleaved: sequences assembled across blocks") {
  const std::string phy =
      "3 6\n"
      "Alpha ACG\n"
      "Beta  TTA\n"
      "Gamma CGC\n"
      "\n"
      "Alpha TAC\n"
      "Beta  GAA\n"
      "Gamma CAT\n";
  std::istringstream ss(phy);
  GeneAlignment gene = read_phylip(ss, "gene", /*interleaved=*/true);

  REQUIRE(gene.num_taxa() == 3);
  REQUIRE(gene.length() == 6);
  REQUIRE(gene.sequences[0].taxon == "Alpha");

  // Alpha: ACG + TAC = ACGTAC
  REQUIRE(gene.sequences[0].bases[0] == NUC_A);
  REQUIRE(gene.sequences[0].bases[1] == NUC_C);
  REQUIRE(gene.sequences[0].bases[2] == NUC_G);
  REQUIRE(gene.sequences[0].bases[3] == NUC_T);
  REQUIRE(gene.sequences[0].bases[4] == NUC_A);
  REQUIRE(gene.sequences[0].bases[5] == NUC_C);
}

TEST_CASE("read_phylip: empty stream throws") {
  std::istringstream ss("");
  REQUIRE_THROWS_AS(read_phylip(ss), std::runtime_error);
}

TEST_CASE("read_phylip: non-numeric header throws") {
  std::istringstream ss("abc\nT1 ACG\n");
  REQUIRE_THROWS_AS(read_phylip(ss), std::runtime_error);
}

// --- read_nexus --------------------------------------------------------------

TEST_CASE("read_nexus: basic sequential DATA block") {
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  DIMENSIONS NTAX=3 NCHAR=4;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  Alpha ACGT\n"
      "  Beta  TTAG\n"
      "  Gamma CGCA\n"
      ";\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  REQUIRE(aln.num_genes() == 1);
  const GeneAlignment& gene = aln.genes[0];
  REQUIRE(gene.num_taxa() == 3);
  REQUIRE(gene.length() == 4);

  auto names = gene.taxon_names();
  REQUIRE(names[0] == "Alpha");
  REQUIRE(names[1] == "Beta");
  REQUIRE(names[2] == "Gamma");

  REQUIRE(gene.sequences[0].bases[0] == NUC_A);
  REQUIRE(gene.sequences[0].bases[3] == NUC_T);
}

TEST_CASE("read_nexus: TITLE statement sets gene_name") {
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  TITLE COI;\n"
      "  DIMENSIONS NTAX=2 NCHAR=3;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  T1 ACG\n"
      "  T2 TTA\n"
      ";\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  REQUIRE(aln.genes[0].gene_name == "COI");
}

TEST_CASE("read_nexus: nested comments stripped, sequence unaffected") {
  const std::string nex =
      "#NEXUS [file-level comment]\n"
      "BEGIN DATA;\n"
      "  DIMENSIONS NTAX=2 NCHAR=3;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  T1 ACG [inline comment]\n"
      "  T2 TTA\n"
      ";\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  REQUIRE(aln.genes[0].sequences[0].bases[0] == NUC_A);
  REQUIRE(aln.genes[0].sequences[0].bases[1] == NUC_C);
  REQUIRE(aln.genes[0].sequences[0].bases[2] == NUC_G);
}

TEST_CASE("read_nexus: multiple DATA blocks yield multiple genes") {
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  TITLE gene_a;\n"
      "  DIMENSIONS NTAX=2 NCHAR=3;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  T1 ACG\n"
      "  T2 TTA\n"
      ";\n"
      "END;\n"
      "BEGIN DATA;\n"
      "  TITLE gene_b;\n"
      "  DIMENSIONS NTAX=2 NCHAR=3;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  T1 GCA\n"
      "  T2 CGT\n"
      ";\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  REQUIRE(aln.num_genes() == 2);
  REQUIRE(aln.genes[0].gene_name == "gene_a");
  REQUIRE(aln.genes[1].gene_name == "gene_b");
  REQUIRE(aln.genes[1].sequences[0].taxon == "T1");
  REQUIRE(aln.genes[1].sequences[0].bases[0] == NUC_G);
}

TEST_CASE("read_nexus: SETS block splits single alignment by charset") {
  // Columns 1-3 → loc1, columns 4-6 → loc2 (1-indexed in NEXUS)
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  DIMENSIONS NTAX=2 NCHAR=6;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  T1 ACGTAC\n"
      "  T2 TTGGCA\n"
      ";\n"
      "END;\n"
      "BEGIN SETS;\n"
      "  CHARSET loc1 = 1-3;\n"
      "  CHARSET loc2 = 4-6;\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  REQUIRE(aln.num_genes() == 2);
  REQUIRE(aln.genes[0].gene_name == "loc1");
  REQUIRE(aln.genes[1].gene_name == "loc2");

  // T1 loc1: positions 0-2 of ACGTAC → A, C, G
  REQUIRE(aln.genes[0].sequences[0].bases[0] == NUC_A);
  REQUIRE(aln.genes[0].sequences[0].bases[1] == NUC_C);
  REQUIRE(aln.genes[0].sequences[0].bases[2] == NUC_G);

  // T1 loc2: positions 3-5 of ACGTAC → T, A, C
  REQUIRE(aln.genes[1].sequences[0].bases[0] == NUC_T);
  REQUIRE(aln.genes[1].sequences[0].bases[1] == NUC_A);
  REQUIRE(aln.genes[1].sequences[0].bases[2] == NUC_C);
}

TEST_CASE("read_nexus: missing #NEXUS header throws") {
  std::istringstream ss("BEGIN DATA;\nEND;\n");
  REQUIRE_THROWS_AS(read_nexus(ss), std::runtime_error);
}

TEST_CASE("read_nexus: unsupported DATATYPE throws") {
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  DIMENSIONS NTAX=2 NCHAR=3;\n"
      "  FORMAT DATATYPE=PROTEIN;\n"
      "  MATRIX\n"
      "  T1 ACG\n"
      "  T2 TTA\n"
      ";\n"
      "END;\n";
  std::istringstream ss(nex);
  REQUIRE_THROWS_AS(read_nexus(ss), std::runtime_error);
}

TEST_CASE("read_nexus: all_taxon_names returns sorted union") {
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  TITLE g1;\n"
      "  DIMENSIONS NTAX=2 NCHAR=2;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  Zebra AC\n"
      "  Ant   TG\n"
      ";\n"
      "END;\n"
      "BEGIN DATA;\n"
      "  TITLE g2;\n"
      "  DIMENSIONS NTAX=2 NCHAR=2;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  Zebra CA\n"
      "  Bear  GT\n"
      ";\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  auto names = aln.all_taxon_names();  // should be sorted
  REQUIRE(names.size() == 3);
  REQUIRE(names[0] == "Ant");
  REQUIRE(names[1] == "Bear");
  REQUIRE(names[2] == "Zebra");
}

// --- compress_sites ----------------------------------------------------------

TEST_CASE("compress_sites: identical columns are collapsed") {
  // 3 taxa, 4 sites: columns 0,2 are (A,A,A); columns 1,3 are (C,G,T)
  GeneAlignment gene = make_gene({{"A", "ACAC"}, {"B", "AGAG"}, {"C", "ATAT"}});
  gene.compress_sites();

  REQUIRE_FALSE(gene.compressed.empty());
  REQUIRE(gene.compressed.n_taxa == 3);
  // (A,A,A) x2, (C,G,T) x2 -> 2 unique patterns
  REQUIRE(gene.compressed.n_patterns == 2);

  // Weights should sum to original length
  double total_weight = 0.0;
  for (double w : gene.compressed.weights) total_weight += w;
  REQUIRE(total_weight == Approx(4.0));
}

TEST_CASE("compress_sites: all identical columns yield one pattern") {
  GeneAlignment gene = make_gene({{"A", "AAAA"}, {"B", "CCCC"}, {"C", "GGGG"}});
  gene.compress_sites();

  REQUIRE(gene.compressed.n_patterns == 1);
  REQUIRE(gene.compressed.weights[0] == Approx(4.0));
}

TEST_CASE("compress_sites: all distinct columns are preserved") {
  GeneAlignment gene = make_gene({{"A", "ACGT"}, {"B", "ACGT"}, {"C", "ACGT"}});
  gene.compress_sites();

  // Each column is unique: (AA,AA,AA), (CC,CC,CC), (GG,GG,GG), (TT,TT,TT)
  REQUIRE(gene.compressed.n_patterns == 4);
  for (size_t i = 0; i < 4; ++i)
    REQUIRE(gene.compressed.weights[i] == Approx(1.0));
}

TEST_CASE("compress_sites: all_defined flags correct with gaps") {
  // Column 0: A,C,G -> all defined
  // Column 1: A,-,G -> has gap
  // Column 2: A,C,G -> all defined (same as col 0, collapsed)
  GeneAlignment gene = make_gene({{"A", "AAA"}, {"B", "C-C"}, {"C", "GGG"}});
  gene.compress_sites();

  REQUIRE(gene.compressed.n_patterns == 2);

  // Find which pattern is the gap one vs the non-gap one
  bool found_gap_pattern = false;
  bool found_nongap_pattern = false;
  for (size_t p = 0; p < gene.compressed.n_patterns; ++p) {
    if (gene.compressed.all_defined[p] == 0) {
      found_gap_pattern = true;
      REQUIRE(gene.compressed.weights[p] == Approx(1.0));
    } else {
      found_nongap_pattern = true;
      REQUIRE(gene.compressed.weights[p] == Approx(2.0));
    }
  }
  REQUIRE(found_gap_pattern);
  REQUIRE(found_nongap_pattern);
}

TEST_CASE("compress_sites: base() accessor returns correct values") {
  GeneAlignment gene = make_gene({{"A", "AC"}, {"B", "GT"}});
  gene.compress_sites();

  REQUIRE(gene.compressed.n_patterns == 2);
  REQUIRE(gene.compressed.n_taxa == 2);

  // Check that we can read the bases back correctly
  // The exact pattern order depends on the map, so check both patterns
  bool found_ag = false, found_ct = false;
  for (size_t p = 0; p < 2; ++p) {
    uint8_t b0 = gene.compressed.base(p, 0);
    uint8_t b1 = gene.compressed.base(p, 1);
    if (b0 == NUC_A && b1 == NUC_G) found_ag = true;
    if (b0 == NUC_C && b1 == NUC_T) found_ct = true;
  }
  REQUIRE(found_ag);
  REQUIRE(found_ct);
}

TEST_CASE("compress_sites: empty alignment stays empty") {
  GeneAlignment gene;
  gene.compress_sites();
  REQUIRE(gene.compressed.empty());
}

TEST_CASE("compress_sites: single-column alignment") {
  GeneAlignment gene = make_gene({{"A", "A"}, {"B", "C"}});
  gene.compress_sites();

  REQUIRE(gene.compressed.n_patterns == 1);
  REQUIRE(gene.compressed.weights[0] == Approx(1.0));
  REQUIRE(gene.compressed.all_defined[0] == 1);
}

// --- compressed vs uncompressed equivalence ----------------------------------

TEST_CASE("compressed and uncompressed paths produce identical site patterns") {
  // Use a non-trivial alignment with duplicate columns and mixed bases
  GeneAlignment gene_raw = make_gene({{"A", "ACGTACGTNN"},
                                       {"B", "TGCATGCA--"},
                                       {"C", "GCATRYNWBH"}});
  // Uncompressed path
  std::vector<TuplePatterns> raw_pats;
  for_each_site_pattern(gene_raw, [&](TuplePatterns p) {
    raw_pats.push_back(std::move(p));
  });
  REQUIRE(raw_pats.size() == 1);

  // Now compress and get site patterns via compressed path
  GeneAlignment gene_comp = make_gene({{"A", "ACGTACGTNN"},
                                        {"B", "TGCATGCA--"},
                                        {"C", "GCATRYNWBH"}});
  gene_comp.compress_sites();
  REQUIRE_FALSE(gene_comp.compressed.empty());

  std::vector<TuplePatterns> comp_pats;
  for_each_site_pattern(gene_comp, [&](TuplePatterns p) {
    comp_pats.push_back(std::move(p));
  });
  REQUIRE(comp_pats.size() == 1);

  // Compare results
  REQUIRE(comp_pats[0].n_total == Approx(raw_pats[0].n_total));
  REQUIRE(comp_pats[0].n_counted == Approx(raw_pats[0].n_counted));
  for (int i = 0; i < 64; ++i) {
    REQUIRE(comp_pats[0].site_pattern[i] ==
            Approx(raw_pats[0].site_pattern[i]).margin(1e-12));
  }
}

TEST_CASE("compressed and uncompressed: 4-taxa ktuple equivalence") {
  GeneAlignment gene_raw = make_gene({{"W", "ACGTACGT"},
                                       {"X", "TGCATGCA"},
                                       {"Y", "GCATGCAT"},
                                       {"Z", "ATCGATCG"}});
  std::vector<TuplePatterns> raw_pats;
  for_each_ktuple_pattern(gene_raw, 4, [&](TuplePatterns p) {
    raw_pats.push_back(std::move(p));
  });
  REQUIRE(raw_pats.size() == 1);

  GeneAlignment gene_comp = make_gene({{"W", "ACGTACGT"},
                                        {"X", "TGCATGCA"},
                                        {"Y", "GCATGCAT"},
                                        {"Z", "ATCGATCG"}});
  gene_comp.compress_sites();

  std::vector<TuplePatterns> comp_pats;
  for_each_ktuple_pattern(gene_comp, 4, [&](TuplePatterns p) {
    comp_pats.push_back(std::move(p));
  });
  REQUIRE(comp_pats.size() == 1);

  REQUIRE(comp_pats[0].n_total == Approx(raw_pats[0].n_total));
  REQUIRE(comp_pats[0].n_counted == Approx(raw_pats[0].n_counted));
  for (size_t i = 0; i < comp_pats[0].site_pattern.size(); ++i) {
    REQUIRE(comp_pats[0].site_pattern[i] ==
            Approx(raw_pats[0].site_pattern[i]).margin(1e-12));
  }
}

TEST_CASE("compressed and uncompressed: site_patterns_for_ktaxa equivalence") {
  GeneAlignment gene_raw = make_gene({{"A", "AACCGGTT"},
                                       {"B", "ACACACAC"},
                                       {"C", "AACCGGTT"}});
  TuplePatterns raw_pat = site_patterns_for_ktaxa(gene_raw, {"A", "B", "C"});

  GeneAlignment gene_comp = make_gene({{"A", "AACCGGTT"},
                                        {"B", "ACACACAC"},
                                        {"C", "AACCGGTT"}});
  gene_comp.compress_sites();
  TuplePatterns comp_pat = site_patterns_for_ktaxa(gene_comp, {"A", "B", "C"});

  REQUIRE(comp_pat.n_total == Approx(raw_pat.n_total));
  REQUIRE(comp_pat.n_counted == Approx(raw_pat.n_counted));
  for (int i = 0; i < 64; ++i) {
    REQUIRE(comp_pat.site_pattern[i] ==
            Approx(raw_pat.site_pattern[i]).margin(1e-12));
  }
}

// --- parsers auto-compress ---------------------------------------------------

TEST_CASE("read_phylip: compressed sites populated after parsing") {
  const std::string phy =
      "3 4\n"
      "Alpha ACGT\n"
      "Beta  TTAG\n"
      "Gamma CGCA\n";
  std::istringstream ss(phy);
  GeneAlignment gene = read_phylip(ss, "mygene");

  REQUIRE_FALSE(gene.compressed.empty());
  REQUIRE(gene.compressed.n_taxa == 3);
  REQUIRE(gene.compressed.n_patterns <= 4);  // at most 4 unique columns
  REQUIRE(gene.compressed.n_patterns >= 1);
}

TEST_CASE("read_nexus: compressed sites populated after parsing") {
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  DIMENSIONS NTAX=3 NCHAR=6;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  T1 AAAAAA\n"
      "  T2 AAAAAA\n"
      "  T3 AAAAAA\n"
      ";\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  REQUIRE(aln.num_genes() == 1);
  REQUIRE_FALSE(aln.genes[0].compressed.empty());
  // All identical columns -> 1 pattern
  REQUIRE(aln.genes[0].compressed.n_patterns == 1);
  REQUIRE(aln.genes[0].compressed.weights[0] == Approx(6.0));
}

TEST_CASE("read_nexus with SETS: each gene partition is compressed") {
  const std::string nex =
      "#NEXUS\n"
      "BEGIN DATA;\n"
      "  DIMENSIONS NTAX=2 NCHAR=6;\n"
      "  FORMAT DATATYPE=DNA;\n"
      "  MATRIX\n"
      "  T1 AAACCC\n"
      "  T2 AAACCC\n"
      ";\n"
      "END;\n"
      "BEGIN SETS;\n"
      "  CHARSET loc1 = 1-3;\n"
      "  CHARSET loc2 = 4-6;\n"
      "END;\n";
  std::istringstream ss(nex);
  MultiGeneAlignment aln = read_nexus(ss);

  REQUIRE(aln.num_genes() == 2);
  // loc1: AAA -> 1 unique pattern with weight 3
  REQUIRE_FALSE(aln.genes[0].compressed.empty());
  REQUIRE(aln.genes[0].compressed.n_patterns == 1);
  REQUIRE(aln.genes[0].compressed.weights[0] == Approx(3.0));
  // loc2: CCC -> 1 unique pattern with weight 3
  REQUIRE_FALSE(aln.genes[1].compressed.empty());
  REQUIRE(aln.genes[1].compressed.n_patterns == 1);
  REQUIRE(aln.genes[1].compressed.weights[0] == Approx(3.0));
}

// --- for_each_site_pattern ---------------------------------------------------

TEST_CASE("for_each_site_pattern: produces n choose 3 triplets") {
  // C(5,3) = 10
  GeneAlignment gene = make_gene({{"A", "ACGT"},
                                  {"B", "ACGT"},
                                  {"C", "ACGT"},
                                  {"D", "ACGT"},
                                  {"E", "ACGT"}});
  int count = 0;
  for_each_site_pattern(gene, [&](TuplePatterns) { ++count; });
  REQUIRE(count == 10);
}

TEST_CASE("for_each_site_pattern: fewer than 3 taxa yields no triplets") {
  for (int n = 0; n <= 2; ++n) {
    std::vector<std::pair<std::string, std::string>> taxa;
    for (int i = 0; i < n; ++i)
      taxa.push_back({std::string(1, static_cast<char>('A' + i)), "ACGT"});
    GeneAlignment gene = make_gene(taxa);
    int count = 0;
    for_each_site_pattern(gene, [&](TuplePatterns) { ++count; });
    REQUIRE(count == 0);
  }
}

TEST_CASE(
    "for_each_site_pattern: taxa in TuplePatterns are lexicographically "
    "sorted") {
  // Input taxa in reverse order: Zebra, Mole, Aardvark
  GeneAlignment gene =
      make_gene({{"Zebra", "ACGT"}, {"Mole", "ACGT"}, {"Aardvark", "ACGT"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.taxa[0] <= pats.taxa[1]);
    REQUIRE(pats.taxa[1] <= pats.taxa[2]);
    REQUIRE(pats.taxa[0] == "Aardvark");
    REQUIRE(pats.taxa[1] == "Mole");
    REQUIRE(pats.taxa[2] == "Zebra");
  });
}

TEST_CASE("for_each_site_pattern: n_total and n_counted") {
  // 4 sites: 3 normal + 1 gap in taxon B
  GeneAlignment gene = make_gene({{"A", "AAAA"}, {"B", "A-AA"}, {"C", "AAAA"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.n_total == Approx(4.0));
    REQUIRE(pats.n_counted == Approx(3.0));
  });
}

TEST_CASE("for_each_site_pattern: all-A sites accumulate in AAA bucket") {
  // Pattern AAA has index 0*16+0*4+0 = 0
  GeneAlignment gene = make_gene({{"A", "AAA"}, {"B", "AAA"}, {"C", "AAA"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.site_pattern[0] == Approx(3.0));
    REQUIRE(total_pattern_mass(pats) == Approx(3.0));
  });
}

TEST_CASE(
    "for_each_site_pattern: distinct unambiguous bases land in correct "
    "bucket") {
  // Sorted taxa: A < B < C
  // A=A(i=0), B=C(j=1), C=G(k=2) → index 0*16+1*4+2 = 6
  GeneAlignment gene = make_gene({{"A", "A"}, {"B", "C"}, {"C", "G"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.n_counted == Approx(1.0));
    REQUIRE(pats.site_pattern[0 * 16 + 1 * 4 + 2] ==
            Approx(1.0));  // A=0, C=1, G=2
    REQUIRE(total_pattern_mass(pats) == Approx(1.0));
  });
}

TEST_CASE(
    "for_each_site_pattern: gap sites are excluded from n and n_counted") {
  // Site 1 has a gap in B; sites 0 and 2 are normal AAA
  GeneAlignment gene = make_gene({{"A", "AAA"}, {"B", "A-A"}, {"C", "AAA"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.n_counted == Approx(2.0));
    REQUIRE(pats.site_pattern[0] == Approx(2.0));  // two AAA sites
    REQUIRE(total_pattern_mass(pats) == Approx(2.0));
  });
}

TEST_CASE(
    "for_each_site_pattern: ambiguous base splits probability fractionally") {
  // Sorted taxa: A < B < C; sequences of length 1
  // A=A(dp={1,0,0,0}), B=C(dq={0,1,0,0}), C=R(dr={0.5,0,0.5,0})
  // Expected: n[ACA]=0.5 (index 0*16+1*4+0=4), n[ACG]=0.5 (index 0*16+1*4+2=6)
  GeneAlignment gene = make_gene({{"A", "A"}, {"B", "C"}, {"C", "R"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.n_counted == Approx(1.0));
    REQUIRE(pats.site_pattern[4] == Approx(0.5));  // ACA
    REQUIRE(pats.site_pattern[6] == Approx(0.5));  // ACG
    REQUIRE(total_pattern_mass(pats) == Approx(1.0));
  });
}

TEST_CASE(
    "for_each_site_pattern: N base distributes equally over all patterns") {
  // A=A, B=A, C=N  →  dp={1,0,0,0}, dq={1,0,0,0}, dr={0.25,0.25,0.25,0.25}
  // AAA=0.25, AAC=0.25, AAG=0.25, AAT=0.25  (indices 0,1,2,3)
  GeneAlignment gene = make_gene({{"A", "A"}, {"B", "A"}, {"C", "N"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.site_pattern[0] == Approx(0.25));  // AAA
    REQUIRE(pats.site_pattern[1] == Approx(0.25));  // AAC
    REQUIRE(pats.site_pattern[2] == Approx(0.25));  // AAG
    REQUIRE(pats.site_pattern[3] == Approx(0.25));  // AAT
    REQUIRE(total_pattern_mass(pats) == Approx(1.0));
  });
}

TEST_CASE(
    "for_each_site_pattern: pattern mass equals n_counted for any sequence") {
  // Invariant: sum of all 64 pattern counts == n_counted (distributions sum to
  // 1)
  GeneAlignment gene = make_gene({{"A", "ACGTRYNWBDHVMSK"},
                                  {"B", "TGCANNRYACGTRYN"},
                                  {"C", "GCATNRYWBDHVMSK"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(total_pattern_mass(pats) == Approx(pats.n_counted));
  });
}

TEST_CASE(
    "for_each_site_pattern: site pattern ordering matches sorted taxon order") {
  // Explicitly verify counts correspond to sorted order, not input order.
  // Input: C="T", A="A", B="C"  (reverse alphabetic)
  // Sorted: A(ip=1)=A, B(ip=2)=C, C(ip=0)=T
  // Site 0: A(i=0), C(j=1), T(k=3) → index 0*16+1*4+3 = 7
  GeneAlignment gene = make_gene({{"C", "T"}, {"A", "A"}, {"B", "C"}});
  for_each_site_pattern(gene, [](TuplePatterns pats) {
    REQUIRE(pats.taxa[0] == "A");
    REQUIRE(pats.taxa[1] == "B");
    REQUIRE(pats.taxa[2] == "C");
    // A=A(0), B=C(1), C=T(3): 0*16 + 1*4 + 3 = 7
    REQUIRE(pats.site_pattern[7] == Approx(1.0));
    REQUIRE(total_pattern_mass(pats) == Approx(1.0));
  });
}

// --- compute_aggregate_site_patterns -----------------------------------------

TEST_CASE("compute_aggregate_site_patterns: accumulates counts across genes") {
  // Gene 1: all-A → AAA (index 0)
  // Gene 2: all-C → CCC (index 1*16+1*4+1 = 21)
  GeneAlignment g1 = make_gene({{"A", "A"}, {"B", "A"}, {"C", "A"}}, "g1");
  GeneAlignment g2 = make_gene({{"A", "C"}, {"B", "C"}, {"C", "C"}}, "g2");

  MultiGeneAlignment aln;
  aln.genes.push_back(g1);
  aln.genes.push_back(g2);

  auto patterns = compute_aggregate_site_patterns(aln);
  REQUIRE(patterns.size() == 1);

  const TuplePatterns& pats = patterns[0];
  REQUIRE(pats.n_total == Approx(2.0));
  REQUIRE(pats.n_counted == Approx(2.0));
  REQUIRE(pats.site_pattern[0] == Approx(1.0));   // AAA from g1
  REQUIRE(pats.site_pattern[21] == Approx(1.0));  // CCC from g2
}

TEST_CASE(
    "compute_aggregate_site_patterns: genes with too few taxa contribute no "
    "triplets") {
  // g1: A,B,C → 1 triplet; g2: A,B only → no triplets from g2
  GeneAlignment g1 = make_gene({{"A", "A"}, {"B", "A"}, {"C", "A"}}, "g1");
  GeneAlignment g2 = make_gene({{"A", "C"}, {"B", "C"}}, "g2");

  MultiGeneAlignment aln;
  aln.genes.push_back(g1);
  aln.genes.push_back(g2);

  auto patterns = compute_aggregate_site_patterns(aln);
  REQUIRE(patterns.size() == 1);
  REQUIRE(patterns[0].n_total == Approx(1.0));  // only g1
  REQUIRE(patterns[0].n_counted == Approx(1.0));
}

TEST_CASE(
    "compute_aggregate_site_patterns: correct triplet count with 4 taxa across "
    "2 genes") {
  // C(4,3) = 4 triplets: (A,B,C), (A,B,D), (A,C,D), (B,C,D)
  GeneAlignment g1 =
      make_gene({{"A", "A"}, {"B", "A"}, {"C", "A"}, {"D", "A"}}, "g1");
  GeneAlignment g2 =
      make_gene({{"A", "C"}, {"B", "C"}, {"C", "C"}, {"D", "C"}}, "g2");

  MultiGeneAlignment aln;
  aln.genes.push_back(g1);
  aln.genes.push_back(g2);

  auto patterns = compute_aggregate_site_patterns(aln);
  REQUIRE(patterns.size() == 4);

  for (const auto& p : patterns) {
    REQUIRE(p.n_counted == Approx(2.0));
    REQUIRE(p.n_total == Approx(2.0));
  }
}
