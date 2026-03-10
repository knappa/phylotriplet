#pragma once
#include <algorithm>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "gen_tree.hpp"  // WeightedTriplet (brings in Eigen3 transitively)
#include "sequence.hpp"
#include "warn.hpp"

namespace phylotriplet {

/** Returns true if the bitmask-encoded nucleotide represents a gap. */
inline bool is_gap(uint8_t b) { return b == NUC_GAP; }

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
  if (b & NUC_A) ++count;
  if (b & NUC_C) ++count;
  if (b & NUC_G) ++count;
  if (b & NUC_T) ++count;
  if (count == 0) return {0.0, 0.0, 0.0, 0.0};
  const double w = 1.0 / count;
  return {(b & NUC_A) ? w : 0.0, (b & NUC_C) ? w : 0.0, (b & NUC_G) ? w : 0.0,
          (b & NUC_T) ? w : 0.0};
}

namespace detail {
/**
 * Recursively accumulate a k-tuple site pattern count.
 *
 * Descends through each of the k positions, multiplying in the base
 * probability at each level and accumulating into site_pattern[idx] at the
 * leaf.  Only recurses into bases with non-zero probability.
 */
inline void accumulate_ktuple(std::vector<double>& site_pattern,
                              const std::vector<NucleotideDistribution>& dists,
                              size_t pos, size_t idx, double weight) {
  if (pos == dists.size()) {
    site_pattern[idx] += weight;
    return;
  }
  for (int b = 0; b < 4; ++b) {
    if (dists[pos][b] > 0.0)
      accumulate_ktuple(site_pattern, dists, pos + 1, idx * 4 + b,
                        weight * dists[pos][b]);
  }
}
}  // namespace detail

/**
 * Call visit(pats) for every triplet in a single gene alignment, producing
 * one TripletSitePatterns at a time without intermediate storage.
 *
 * Taxa within each TripletSitePatterns are in lexicographic order so that
 * results from different genes can be aggregated by key.
 *
 * site_pattern[i*16 + j*4 + k] counts (fractionally) sites where taxa[0]
 * shows base i, taxa[1] base j, and taxa[2] base k (0=A,1=C,2=G,3=T).
 * Ambiguous bases contribute fractionally to all compatible patterns.
 */
template <typename Visitor>
void for_each_site_pattern(const GeneAlignment& gene, Visitor&& visit) {
  const size_t n = gene.num_taxa();
  const size_t len = gene.length();

  if (n < 3) return;

  // Build a sorted index so taxa come out in lexicographic order.
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return gene.sequences[a].taxon < gene.sequences[b].taxon;
  });

  for (size_t p = 0; p < n; ++p) {
    for (size_t q = p + 1; q < n; ++q) {
      for (size_t r = q + 1; r < n; ++r) {
        const size_t ip = order[p], iq = order[q], ir = order[r];

        TuplePatterns pats({gene.sequences[ip].taxon, gene.sequences[iq].taxon,
                            gene.sequences[ir].taxon});
        pats.n_total = static_cast<double>(len);

        const std::vector<uint8_t>& sp = gene.sequences[ip].bases;
        const std::vector<uint8_t>& sq = gene.sequences[iq].bases;
        const std::vector<uint8_t>& sr = gene.sequences[ir].bases;

        for (size_t s = 0; s < len; ++s) {
          const uint8_t bp = (s < sp.size()) ? sp[s] : NUC_GAP;
          const uint8_t bq = (s < sq.size()) ? sq[s] : NUC_GAP;
          const uint8_t br = (s < sr.size()) ? sr[s] : NUC_GAP;

          if (is_gap(bp) || is_gap(bq) || is_gap(br)) continue;

          pats.n_counted += 1.0;

          const NucleotideDistribution dp = nucleotide_distribution(bp);
          const NucleotideDistribution dq = nucleotide_distribution(bq);
          const NucleotideDistribution dr = nucleotide_distribution(br);

          for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
              for (int k = 0; k < 4; ++k)
                pats.site_pattern[i * 16 + j * 4 + k] += dp[i] * dq[j] * dr[k];
        }

        visit(std::move(pats));
      }
    }
  }
}

/**
 * Call visit(pats) for every k-tuple of taxa in a gene alignment, producing
 * one TuplePatterns at a time without intermediate storage.  k must be >= 2.
 *
 * Taxa within each TuplePatterns are in lexicographic order.
 * site_pattern has 4^k entries indexed as described in TuplePatterns.
 */
template <typename Visitor>
void for_each_ktuple_pattern(const GeneAlignment& gene, size_t k,
                             Visitor&& visit) {
  if (k < 2)
    throw std::invalid_argument("for_each_ktuple_pattern: k must be >= 2");

  const size_t n = gene.num_taxa();
  const size_t len = gene.length();
  if (n < k) return;

  // Build a sorted index so taxa come out in lexicographic order.
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return gene.sequences[a].taxon < gene.sequences[b].taxon;
  });

  // Iterate over all C(n, k) combinations using an index vector.
  std::vector<size_t> combo(k);
  for (size_t i = 0; i < k; ++i) combo[i] = i;

  while (true) {
    std::vector<std::string> taxa(k);
    std::vector<const std::vector<uint8_t>*> seqs(k);
    for (size_t t = 0; t < k; ++t) {
      const size_t idx = order[combo[t]];
      taxa[t] = gene.sequences[idx].taxon;
      seqs[t] = &gene.sequences[idx].bases;
    }

    TuplePatterns pats(std::move(taxa));
    pats.n_total = static_cast<double>(len);

    for (size_t s = 0; s < len; ++s) {
      std::vector<NucleotideDistribution> dists(k);
      bool any_gap = false;
      for (size_t t = 0; t < k; ++t) {
        const uint8_t b = (s < seqs[t]->size()) ? (*seqs[t])[s] : NUC_GAP;
        if (is_gap(b)) {
          any_gap = true;
          break;
        }
        dists[t] = nucleotide_distribution(b);
      }
      if (any_gap) continue;
      pats.n_counted += 1.0;
      detail::accumulate_ktuple(pats.site_pattern, dists, 0, 0, 1.0);
    }

    visit(std::move(pats));

    // Advance to next combination in lexicographic order.
    int pos = static_cast<int>(k) - 1;
    while (pos >= 0 && combo[pos] == n - k + static_cast<size_t>(pos)) --pos;
    if (pos < 0) break;
    ++combo[pos];
    for (size_t i = static_cast<size_t>(pos) + 1; i < k; ++i)
      combo[i] = combo[i - 1] + 1;
  }
}

/**
 * Compute site-pattern counts for a specific k-tuple of taxa in one gene.
 *
 * taxa may be in any order; they are sorted lexicographically internally.
 * k must be >= 2.  Throws std::invalid_argument if any taxon is missing or
 * names are duplicated.
 */
inline TuplePatterns site_patterns_for_ktaxa(const GeneAlignment& gene,
                                             std::vector<std::string> taxa) {
  std::sort(taxa.begin(), taxa.end());
  if (taxa.size() < 2)
    throw std::invalid_argument(
        "site_patterns_for_ktaxa: need at least 2 taxa");
  for (size_t t = 0; t + 1 < taxa.size(); ++t)
    if (taxa[t] == taxa[t + 1])
      throw std::invalid_argument(
          "site_patterns_for_ktaxa: duplicate taxon names");

  const size_t k = taxa.size();
  std::vector<const std::vector<uint8_t>*> seqs(k, nullptr);
  for (const auto& seq : gene.sequences)
    for (size_t t = 0; t < k; ++t)
      if (seq.taxon == taxa[t]) seqs[t] = &seq.bases;

  for (size_t t = 0; t < k; ++t)
    if (seqs[t] == nullptr)
      throw std::invalid_argument("site_patterns_for_ktaxa: taxon '" + taxa[t] +
                                  "' not found in alignment '" +
                                  gene.gene_name + "'");

  TuplePatterns pats(taxa);
  const size_t len = gene.length();
  pats.n_total = static_cast<double>(len);

  for (size_t s = 0; s < len; ++s) {
    std::vector<NucleotideDistribution> dists(k);
    bool any_gap = false;
    for (size_t t = 0; t < k; ++t) {
      const uint8_t b = (s < seqs[t]->size()) ? (*seqs[t])[s] : NUC_GAP;
      if (is_gap(b)) {
        any_gap = true;
        break;
      }
      dists[t] = nucleotide_distribution(b);
    }
    if (any_gap) continue;
    pats.n_counted += 1.0;
    detail::accumulate_ktuple(pats.site_pattern, dists, 0, 0, 1.0);
  }

  return pats;
}

/**
 * Compute aggregate site-pattern counts for a specific k-tuple of taxa
 * across all genes.  
 * k must be >= 2.  Throws if names are duplicated or any taxon is absent
 * from every gene.
 */
inline TuplePatterns aggregate_site_patterns_for_ktaxa(
    const MultiGeneAlignment& alignment, std::vector<std::string> taxa) {
  std::sort(taxa.begin(), taxa.end());
  if (taxa.size() < 2)
    throw std::invalid_argument(
        "aggregate_site_patterns_for_ktaxa: need at least 2 taxa");
  for (size_t t = 0; t + 1 < taxa.size(); ++t)
    if (taxa[t] == taxa[t + 1])
      throw std::invalid_argument(
          "aggregate_site_patterns_for_ktaxa: duplicate taxon names");

  const size_t k = taxa.size();

  for (size_t t = 0; t < k; ++t) {
    bool found = false;
    for (const auto& gene : alignment.genes)
      for (const auto& seq : gene.sequences)
        if (seq.taxon == taxa[t]) {
          found = true;
          break;
        }
    if (!found)
      throw std::invalid_argument("aggregate_site_patterns_for_ktaxa: taxon '" +
                                  taxa[t] + "' not found in any gene");
  }

  TuplePatterns result(taxa);
  int n_skipped = 0;

  for (const auto& gene : alignment.genes) {
    bool all_present = true;
    for (size_t t = 0; t < k; ++t) {
      bool found = false;
      for (const auto& seq : gene.sequences)
        if (seq.taxon == taxa[t]) {
          found = true;
          break;
        }
      if (!found) {
        all_present = false;
        break;
      }
    }
    if (!all_present) {
      ++n_skipped;
      continue;
    }

    TuplePatterns pats = site_patterns_for_ktaxa(gene, taxa);
    for (size_t i = 0; i < result.site_pattern.size(); ++i)
      result.site_pattern[i] += pats.site_pattern[i];
    result.n_total += pats.n_total;
    result.n_counted += pats.n_counted;
  }

  if (n_skipped > 0)
    PHYLOTRIPLET_WARN(std::string("aggregate_site_patterns: ") +
                      std::to_string(n_skipped) + " of " +
                      std::to_string(alignment.num_genes()) +
                      " gene(s) were skipped (did not contain all " +
                      std::to_string(k) + " requested taxa)");

  return result;
}

/**
 * Compute site-pattern counts for a specific taxon triplet within one gene.
 * Convenience wrapper around site_patterns_for_ktaxa.
 */
inline TripletSitePatterns site_patterns_for_taxa(const GeneAlignment& gene,
                                                  const std::string& taxon_a,
                                                  const std::string& taxon_b,
                                                  const std::string& taxon_c) {
  return site_patterns_for_ktaxa(gene, {taxon_a, taxon_b, taxon_c});
}

/**
 * Convenience wrapper: collect all TripletSitePatterns for a gene into a
 * vector. The result has n choose 3 elements; for large n prefer
 * for_each_site_pattern.
 */
inline std::vector<TripletSitePatterns> compute_site_patterns(
    const GeneAlignment& gene) {
  std::vector<TripletSitePatterns> result;
  const size_t n = gene.num_taxa();
  if (n >= 3) result.reserve(n * (n - 1) * (n - 2) / 6);
  for_each_site_pattern(gene, [&](TripletSitePatterns pats) {
    result.push_back(std::move(pats));
  });
  return result;
}

/**
 * Compute aggregate site pattern counts by streaming over all genes.
 *
 * For each triplet, contributions are included only from genes that contain
 * all three taxa.  n_total and n_counted are summed accordingly.
 * Only one gene's worth of per-triplet data is live at a time.
 */
inline std::vector<TripletSitePatterns> compute_aggregate_site_patterns(
    const MultiGeneAlignment& alignment) {
  using Key = std::tuple<std::string, std::string, std::string>;
  std::map<Key, TripletSitePatterns> acc;
  int n_skipped = 0;

  for (const auto& gene : alignment.genes) {
    if (gene.num_taxa() < 3) {
      ++n_skipped;
      continue;
    }
    for_each_site_pattern(gene, [&](TripletSitePatterns pats) {
      Key key(pats.taxa[0], pats.taxa[1], pats.taxa[2]);
      TripletSitePatterns& a = acc[key];
      if (a.taxa.empty()) {
        a.taxa = pats.taxa;
        a.site_pattern.assign(pats.site_pattern.size(), 0.0);
      }
      for (size_t i = 0; i < a.site_pattern.size(); ++i)
        a.site_pattern[i] += pats.site_pattern[i];
      a.n_total += pats.n_total;
      a.n_counted += pats.n_counted;
    });
  }

  if (n_skipped > 0)
    PHYLOTRIPLET_WARN(std::string("triplet_site_patterns: ") +
                      std::to_string(n_skipped) + " of " +
                      std::to_string(alignment.num_genes()) +
                      " gene(s) were skipped (fewer than 3 taxa)");

  std::vector<TripletSitePatterns> result;
  result.reserve(acc.size());
  for (auto& kv : acc) result.push_back(std::move(kv.second));
  return result;
}

/**
 * Compute aggregate site-pattern counts for a specific taxon triplet across
 * all genes in a MultiGeneAlignment.
 *
 * The three names may be given in any order; they are sorted lexicographically
 * internally.
 *
 * Throws std::invalid_argument if:
 *   - the same taxon name is given more than once, or
 *   - any of the three taxa is absent from every gene in the alignment.
 */
/**
 * Compute aggregate site-pattern counts for a specific taxon triplet.
 * Convenience wrapper around aggregate_site_patterns_for_ktaxa.
 */
inline TripletSitePatterns aggregate_site_patterns_for_taxa(
    const MultiGeneAlignment& alignment, const std::string& taxon_a,
    const std::string& taxon_b, const std::string& taxon_c) {
  return aggregate_site_patterns_for_ktaxa(alignment,
                                           {taxon_a, taxon_b, taxon_c});
}

/**
 * STUB: Compute weights for the three possible rooted triplet topologies for 
 * a specific taxon triplet, based on site-pattern counts across all genes.
 *
 * Not yet implemented — always throws std::runtime_error.
 */
inline std::array<double, 3> compute_triplet_weights(
    const MultiGeneAlignment& alignment, const std::string& taxon_a,
    const std::string& taxon_b, const std::string& taxon_c) {
  throw std::runtime_error(
      "compute_triplet_weights: not yet implemented");
}

/**
 * STUB: Compute weights for the three possible rooted triplet topologies for 
 * a specific taxon triplet, based on site-pattern counts in a single gene.
 *
 * Not yet implemented — always throws std::runtime_error.
 */
inline std::array<double, 3> compute_triplet_weights(
    const GeneAlignment& alignment, const std::string& taxon_a,
    const std::string& taxon_b, const std::string& taxon_c) {
  throw std::runtime_error(
      "compute_triplet_weights: not yet implemented");
}


// ---------------------------------------------------------------------------
// Named public API
// ---------------------------------------------------------------------------

/**
 * Compute aggregate site-pattern counts for all C(n,3) taxon triplets across
 * all genes.  Equivalent to compute_aggregate_site_patterns().
 */
inline std::vector<TripletSitePatterns> triplet_site_patterns(
    const MultiGeneAlignment& alignment) {
  return compute_aggregate_site_patterns(alignment);
}

/**
 * Compute aggregate site-pattern counts for a specific taxon triplet across
 * all genes.  Equivalent to aggregate_site_patterns_for_taxa().
 */
inline TripletSitePatterns triplet_site_pattern(
    const MultiGeneAlignment& alignment, const std::string& taxon_a,
    const std::string& taxon_b, const std::string& taxon_c) {
  return aggregate_site_patterns_for_taxa(alignment, taxon_a, taxon_b, taxon_c);
}

/**
 * Compute aggregate site-pattern counts for all C(n,4) taxon quartets across
 * all genes.  Only genes containing all four taxa contribute.
 */
inline std::vector<TuplePatterns> quartet_site_patterns(
    const MultiGeneAlignment& alignment) {
  using Key = std::tuple<std::string, std::string, std::string, std::string>;
  std::map<Key, TuplePatterns> acc;
  int n_skipped = 0;

  for (const auto& gene : alignment.genes) {
    if (gene.num_taxa() < 4) {
      ++n_skipped;
      continue;
    }
    for_each_ktuple_pattern(gene, 4, [&](TuplePatterns pats) {
      Key key(pats.taxa[0], pats.taxa[1], pats.taxa[2], pats.taxa[3]);
      TuplePatterns& a = acc[key];
      if (a.taxa.empty()) {
        a = TuplePatterns(pats.taxa);
      }
      for (size_t i = 0; i < a.site_pattern.size(); ++i)
        a.site_pattern[i] += pats.site_pattern[i];
      a.n_total += pats.n_total;
      a.n_counted += pats.n_counted;
    });
  }

  if (n_skipped > 0)
    PHYLOTRIPLET_WARN(std::string("quartet_site_patterns: ") +
                      std::to_string(n_skipped) + " of " +
                      std::to_string(alignment.num_genes()) +
                      " gene(s) were skipped (fewer than 4 taxa)");

  std::vector<TuplePatterns> result;
  result.reserve(acc.size());
  for (auto& kv : acc) result.push_back(std::move(kv.second));
  return result;
}

/**
 * Compute aggregate site-pattern counts for a specific taxon quartet across
 * all genes.  Equivalent to aggregate_site_patterns_for_ktaxa() with k=4.
 */
inline TuplePatterns quartet_site_pattern(const MultiGeneAlignment& alignment,
                                          const std::string& taxon_a,
                                          const std::string& taxon_b,
                                          const std::string& taxon_c,
                                          const std::string& taxon_d) {
  return aggregate_site_patterns_for_ktaxa(
      alignment, {taxon_a, taxon_b, taxon_c, taxon_d});
}

/**
 * Compute aggregate site-pattern counts for a user-specified k-tuple of taxa
 * (k >= 2) across all genes.  Equivalent to
 * aggregate_site_patterns_for_ktaxa().
 */
inline TuplePatterns site_pattern(const MultiGeneAlignment& alignment,
                                  std::vector<std::string> taxa) {
  return aggregate_site_patterns_for_ktaxa(alignment, std::move(taxa));
}

// ---------------------------------------------------------------------------
// Single-gene overloads (no aggregation, no skipping)
// ---------------------------------------------------------------------------

/** All C(n,3) triplet site-pattern counts for one gene. */
inline std::vector<TripletSitePatterns> triplet_site_patterns(
    const GeneAlignment& gene) {
  return compute_site_patterns(gene);
}

/** Site-pattern counts for a specific triplet within one gene.
 *  Throws std::invalid_argument if any taxon is absent. */
inline TripletSitePatterns triplet_site_pattern(const GeneAlignment& gene,
                                                const std::string& taxon_a,
                                                const std::string& taxon_b,
                                                const std::string& taxon_c) {
  return site_patterns_for_taxa(gene, taxon_a, taxon_b, taxon_c);
}

/** All C(n,4) quartet site-pattern counts for one gene. */
inline std::vector<TuplePatterns> quartet_site_patterns(
    const GeneAlignment& gene) {
  std::vector<TuplePatterns> result;
  for_each_ktuple_pattern(
      gene, 4, [&](TuplePatterns p) { result.push_back(std::move(p)); });
  return result;
}

/** Site-pattern counts for a specific quartet within one gene.
 *  Throws std::invalid_argument if any taxon is absent. */
inline TuplePatterns quartet_site_pattern(const GeneAlignment& gene,
                                          const std::string& taxon_a,
                                          const std::string& taxon_b,
                                          const std::string& taxon_c,
                                          const std::string& taxon_d) {
  return site_patterns_for_ktaxa(gene, {taxon_a, taxon_b, taxon_c, taxon_d});
}

/** Site-pattern counts for a specific k-tuple within one gene.
 *  Throws std::invalid_argument if any taxon is absent. */
inline TuplePatterns site_pattern(const GeneAlignment& gene,
                                  std::vector<std::string> taxa) {
  return site_patterns_for_ktaxa(gene, std::move(taxa));
}

}  // namespace phylotriplet