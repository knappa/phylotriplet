#pragma once
#include <array>
#include <cctype>
#include <cstdint>
#include <map>
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

/** Returns true if the bitmask-encoded nucleotide represents a gap. */
inline bool is_gap(uint8_t b) {
  return b == NUC_GAP;
}

/**
 * Map an unambiguous bitmask to its base index (A=0, C=1, G=2, T=3).
 * Returns -1 for gap or ambiguous bases.
 */
inline int unambiguous_base_index(uint8_t b) {
  if (b == NUC_A)
    return 0;
  if (b == NUC_C)
    return 1;
  if (b == NUC_G)
    return 2;
  if (b == NUC_T)
    return 3;
  return -1;
}

/**
 * Map a bitmask-encoded nucleotide to its probability distribution over
 * {A=0, C=1, G=2, T=3}.
 *
 * Ambiguity codes are treated as equiprobable mixtures over the bases they
 * represent (e.g. R = 0.5 A + 0.5 G).  NUC_GAP returns {0,0,0,0}; callers
 * should check is_gap() first.
 */
inline NucleotideDistribution nucleotide_distribution(uint8_t b) {
  int count = 0;
  if (b & NUC_A)
    ++count;
  if (b & NUC_C)
    ++count;
  if (b & NUC_G)
    ++count;
  if (b & NUC_T)
    ++count;
  if (count == 0)
    return {0.0, 0.0, 0.0, 0.0};
  const double w = 1.0 / count;
  return {
      (b & NUC_A) ? w : 0.0, (b & NUC_C) ? w : 0.0, (b & NUC_G) ? w : 0.0, (b & NUC_T) ? w : 0.0};
}

/** A single named sequence; bases are bitmask-encoded (see NUC_* constants). */
struct Sequence {
  std::string taxon;
  std::vector<uint8_t> bases;
};

/**
 * Column-deduplicated site data for a GeneAlignment.
 *
 * Each unique column pattern (across all taxa) is stored once with a weight
 * equal to the number of times it appears.  Data is column-major:
 * pattern_data[p * n_taxa + t] is the base for taxon t in pattern p.
 */
struct CompressedSites {
  std::vector<uint8_t> pattern_data;  // P * n_taxa entries, column-major
  std::vector<double> weights;        // P pattern counts
  std::vector<uint8_t> all_defined;   // P flags: 1 = no taxon has a gap
  size_t n_patterns = 0;
  size_t n_taxa = 0;

  bool empty() const { return n_patterns == 0; }

  uint8_t base(size_t pattern, size_t taxon) const {
    return pattern_data[pattern * n_taxa + taxon];
  }
};

/** A single-gene alignment: all sequences must have the same length. */
struct GeneAlignment {
  std::string gene_name;
  std::vector<Sequence> sequences;

  size_t num_taxa() const { return sequences.size(); }
  size_t length() const { return sequences.empty() ? 0 : sequences[0].bases.size(); }

  std::vector<std::string> taxon_names() const {
    std::vector<std::string> ns;
    ns.reserve(sequences.size());
    for (const auto& s : sequences)
      ns.push_back(s.taxon);
    return ns;
  }

  /**
   * Empirical nucleotide frequencies (pi_A, pi_C, pi_G, pi_T).
   * Gaps are excluded; ambiguous bases contribute fractionally.
   * Returns {0,0,0,0} for an all-gap or empty alignment.
   */
  NucleotideDistribution nucleotide_frequencies() const {
    NucleotideDistribution freq = {0.0, 0.0, 0.0, 0.0};
    double total = 0.0;
    for (const auto& seq : sequences) {
      for (uint8_t b : seq.bases) {
        if (is_gap(b))
          continue;
        const NucleotideDistribution d = nucleotide_distribution(b);
        for (int i = 0; i < 4; ++i)
          freq[i] += d[i];
        total += 1.0;
      }
    }
    if (total > 0.0)
      for (int i = 0; i < 4; ++i)
        freq[i] /= total;
    return freq;
  }

  /**
   * Deduplicated column data.  Empty until compress_sites() is called.
   * Site pattern functions use this automatically when available.
   */
  CompressedSites compressed;

  /**
   * Build a column-deduplicated representation of the alignment.
   *
   * Identical columns (across all taxa) are collapsed into unique patterns
   * with weights.  Each pattern also records whether it is gap-free for all
   * taxa.  Call after filter_invariant_sites() for maximum benefit.
   *
   * Safe to call multiple times (rebuilds from current sequence data).
   */
  void compress_sites() {
    if (sequences.empty())
      return;
    const size_t n = num_taxa();
    const size_t len = length();
    if (len == 0)
      return;

    compressed = CompressedSites();
    compressed.n_taxa = n;

    std::map<std::vector<uint8_t>, size_t> col_map;
    std::vector<uint8_t> col(n);

    for (size_t s = 0; s < len; ++s) {
      for (size_t t = 0; t < n; ++t) {
        col[t] = sequences[t].bases[s];
      }

      auto it = col_map.find(col);
      if (it != col_map.end()) {
        compressed.weights[it->second] += 1.0;
      } else {
        const size_t idx = compressed.n_patterns++;
        col_map[col] = idx;
        compressed.pattern_data.insert(compressed.pattern_data.end(), col.begin(), col.end());
        compressed.weights.push_back(1.0);

        uint8_t defined = 1;
        for (size_t t = 0; t < n; ++t) {
          if (is_gap(col[t])) {
            defined = 0;
            break;
          }
        }
        compressed.all_defined.push_back(defined);
      }
    }
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
    for (const auto& g : genes)
      ns.push_back(g.gene_name);
    return ns;
  }

  /** Sorted union of all taxon names across all genes. */
  std::vector<std::string> all_taxon_names() const {
    std::set<std::string> seen;
    for (const auto& g : genes)
      for (const auto& s : g.sequences)
        seen.insert(s.taxon);
    return {seen.begin(), seen.end()};
  }

  /**
   * Empirical nucleotide frequencies (pi_A, pi_C, pi_G, pi_T) across all
   * genes.  Gaps are excluded; ambiguous bases contribute fractionally.
   * Returns {0,0,0,0} if all sites are gaps or the alignment is empty.
   */
  NucleotideDistribution nucleotide_frequencies() const {
    NucleotideDistribution freq = {0.0, 0.0, 0.0, 0.0};
    double total = 0.0;
    for (const auto& gene : genes) {
      for (const auto& seq : gene.sequences) {
        for (uint8_t b : seq.bases) {
          if (is_gap(b))
            continue;
          const NucleotideDistribution d = nucleotide_distribution(b);
          for (int i = 0; i < 4; ++i)
            freq[i] += d[i];
          total += 1.0;
        }
      }
    }
    if (total > 0.0)
      for (int i = 0; i < 4; ++i)
        freq[i] /= total;
    return freq;
  }
};

namespace detail {
/// Compute 4^k at runtime.
inline size_t pow4_rt(size_t k) {
  size_t r = 1;
  for (size_t i = 0; i < k; ++i)
    r *= 4;
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
      : taxa(std::move(sorted_taxa)), site_pattern(detail::pow4_rt(taxa.size()), 0.0) {}

  /** Number of taxa in the tuple. */
  size_t order() const { return taxa.size(); }
};

}  // namespace phylotriplet
