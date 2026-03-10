#pragma once
#include <array>
#include <cctype>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

namespace phylotriplet {

/**
 * Nucleotide probability distribution over {A=0, C=1, G=2, T=3}.
 * U is mapped to T for RNA/DNA equivalence.
 */
using NucleotideDistribution = std::array<double, 4>;

/**
 * Compact nucleotide bitmask encoding.
 * Each nucleotide occupies one bit; ambiguity codes are OR-combinations.
 * Gap / unknown = 0x0.
 *
 *   Bit 0 (0x1) = A
 *   Bit 1 (0x2) = C
 *   Bit 2 (0x4) = G
 *   Bit 3 (0x8) = T  (U is also mapped to T)
 */
constexpr uint8_t NUC_A = 0x1;
constexpr uint8_t NUC_C = 0x2;
constexpr uint8_t NUC_G = 0x4;
constexpr uint8_t NUC_T = 0x8;
constexpr uint8_t NUC_GAP = 0x0;

// IUPAC ambiguity codes as bitmask combinations
constexpr uint8_t NUC_R = NUC_A | NUC_G;                  // A or G (purine)
constexpr uint8_t NUC_Y = NUC_C | NUC_T;                  // C or T (pyrimidine)
constexpr uint8_t NUC_S = NUC_C | NUC_G;                  // C or G (strong)
constexpr uint8_t NUC_W = NUC_A | NUC_T;                  // A or T (weak)
constexpr uint8_t NUC_K = NUC_G | NUC_T;                  // G or T (keto)
constexpr uint8_t NUC_M = NUC_A | NUC_C;                  // A or C (amino)
constexpr uint8_t NUC_B = NUC_C | NUC_G | NUC_T;          // not A
constexpr uint8_t NUC_D = NUC_A | NUC_G | NUC_T;          // not C
constexpr uint8_t NUC_H = NUC_A | NUC_C | NUC_T;          // not G
constexpr uint8_t NUC_V = NUC_A | NUC_C | NUC_G;          // not T
constexpr uint8_t NUC_N = NUC_A | NUC_C | NUC_G | NUC_T;  // any

/**
 * Convert an IUPAC nucleotide character (case-insensitive) to a bitmask.
 * U is mapped to T.  Gap characters ('-', '?') and unrecognised characters
 * return NUC_GAP (0).
 */
inline uint8_t char_to_nucleotide(char c) {
  switch (std::toupper(static_cast<unsigned char>(c))) {
    case 'A':
      return NUC_A;
    case 'C':
      return NUC_C;
    case 'G':
      return NUC_G;
    case 'T':
    case 'U':
      return NUC_T;
    case 'R':
      return NUC_R;
    case 'Y':
      return NUC_Y;
    case 'S':
      return NUC_S;
    case 'W':
      return NUC_W;
    case 'K':
      return NUC_K;
    case 'M':
      return NUC_M;
    case 'B':
      return NUC_B;
    case 'D':
      return NUC_D;
    case 'H':
      return NUC_H;
    case 'V':
      return NUC_V;
    case 'N':
    case 'X':
      return NUC_N;
    default:
      return NUC_GAP;
  }
}

/** A single named sequence; bases are bitmask-encoded (see NUC_* constants). */
struct Sequence {
  std::string taxon;
  std::vector<uint8_t> bases;
};

/** A single-gene alignment: all sequences must have the same length. */
struct GeneAlignment {
  std::string gene_name;
  std::vector<Sequence> sequences;

  size_t num_taxa() const { return sequences.size(); }
  size_t length() const {
    return sequences.empty() ? 0 : sequences[0].bases.size();
  }

  std::vector<std::string> taxon_names() const {
    std::vector<std::string> ns;
    ns.reserve(sequences.size());
    for (const auto& s : sequences) ns.push_back(s.taxon);
    return ns;
  }
};

/**
 * A collection of gene alignments.
 * Taxa may differ between genes (partial data is allowed).
 */
struct MultiGeneAlignment {
  std::vector<GeneAlignment> genes;

  size_t num_genes() const { return genes.size(); }

  /** Gene names in order. */
  std::vector<std::string> gene_names() const {
    std::vector<std::string> ns;
    ns.reserve(genes.size());
    for (const auto& g : genes) ns.push_back(g.gene_name);
    return ns;
  }

  /** Sorted union of all taxon names across all genes. */
  std::vector<std::string> all_taxon_names() const {
    std::set<std::string> seen;
    for (const auto& g : genes)
      for (const auto& s : g.sequences) seen.insert(s.taxon);
    return {seen.begin(), seen.end()};
  }
};

namespace detail {
/// Compute 4^k at runtime.
inline size_t pow4_rt(size_t k) {
  size_t r = 1;
  for (size_t i = 0; i < k; ++i) r *= 4;
  return r;
}
}  // namespace detail

/**
 * Site pattern counts for an ordered k-tuple of taxa (k >= 2).
 *
 * Taxa are stored in lexicographic order so that the same biological tuple
 * always maps to the same struct regardless of which gene it came from.
 *
 * site_pattern has 4^k entries.  The flat index for a column where taxa[0]
 * shows base b0, taxa[1] shows b1, ..., taxa[k-1] shows b_{k-1} is:
 *   b0 * 4^(k-1) + b1 * 4^(k-2) + ... + b_{k-1}
 * with 0=A, 1=C, 2=G, 3=T.  Ambiguous bases contribute fractionally to all
 * compatible patterns (e.g. R = 0.5*A + 0.5*G).
 *
 * Gap columns (any taxon has NUC_GAP) are excluded from site_pattern and
 * n_counted but are included in n_total.
 */
struct TuplePatterns {
  std::vector<std::string> taxa;     // sorted lexicographically; length = k
  std::vector<double> site_pattern;  // 4^k entries, zero-initialised
  double n_total = 0.0;
  double n_counted = 0.0;

  TuplePatterns() = default;

  /** Construct for k taxa that are already sorted lexicographically. */
  explicit TuplePatterns(std::vector<std::string> sorted_taxa)
      : taxa(std::move(sorted_taxa)),
        site_pattern(detail::pow4_rt(taxa.size()), 0.0) {}

  /** Number of taxa in the tuple. */
  size_t order() const { return taxa.size(); }
};

/** Backward-compatible alias: a 3-taxon TuplePatterns. */
using TripletSitePatterns = TuplePatterns;

}  // namespace phylotriplet
