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
#include "triplet_weights.hpp"
#include "warn.hpp"

namespace phylotriplet {

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
                              size_t pos,
                              size_t idx,
                              double weight) {
  if (pos == dists.size()) {
    site_pattern[idx] += weight;
    return;
  }
  for (int b = 0; b < 4; ++b) {
    if (dists[pos][b] > 0.0)
      accumulate_ktuple(
          site_pattern, dists, pos + 1, idx * 4 + b, weight * dists[pos][b]);
  }
}

}  // namespace detail

/**
 * Call visit(pats) for every triplet in a single gene alignment, producing
 * one TuplePatterns at a time without intermediate storage.
 *
 * Taxa within each TuplePatterns are in lexicographic order so that
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

  if (n < 3)
    return;

  // Build a sorted index so taxa come out in lexicographic order.
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return gene.sequences[a].taxon < gene.sequences[b].taxon;
  });

  const bool use_compressed = !gene.compressed.empty();

  for (size_t p = 0; p < n; ++p) {
    for (size_t q = p + 1; q < n; ++q) {
      for (size_t r = q + 1; r < n; ++r) {
        const size_t ip = order[p], iq = order[q], ir = order[r];

        TuplePatterns pats({gene.sequences[ip].taxon,
                            gene.sequences[iq].taxon,
                            gene.sequences[ir].taxon});
        pats.n_total = static_cast<double>(gene.length());

        if (use_compressed) {
          const auto& cs = gene.compressed;
          for (size_t u = 0; u < cs.n_patterns; ++u) {
            const uint8_t bp = cs.base(u, ip);
            const uint8_t bq = cs.base(u, iq);
            const uint8_t br = cs.base(u, ir);

            if (!cs.all_defined[u] && (is_gap(bp) || is_gap(bq) || is_gap(br)))
              continue;

            const double w = cs.weights[u];
            pats.n_counted += w;

            const int i2 = unambiguous_base_index(bp);
            const int j2 = unambiguous_base_index(bq);
            const int k2 = unambiguous_base_index(br);
            if (i2 >= 0 && j2 >= 0 && k2 >= 0) {
              pats.site_pattern[i2 * 16 + j2 * 4 + k2] += w;
            } else {
              const NucleotideDistribution dp = nucleotide_distribution(bp);
              const NucleotideDistribution dq = nucleotide_distribution(bq);
              const NucleotideDistribution dr = nucleotide_distribution(br);
              for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                  for (int k = 0; k < 4; ++k)
                    pats.site_pattern[i * 16 + j * 4 + k] +=
                        dp[i] * dq[j] * dr[k] * w;
            }
          }
        } else {
          const std::vector<uint8_t>& sp = gene.sequences[ip].bases;
          const std::vector<uint8_t>& sq = gene.sequences[iq].bases;
          const std::vector<uint8_t>& sr = gene.sequences[ir].bases;

          for (size_t s = 0; s < len; ++s) {
            const uint8_t bp = sp[s], bq = sq[s], br = sr[s];

            if (is_gap(bp) || is_gap(bq) || is_gap(br))
              continue;

            pats.n_counted += 1.0;

            const int i2 = unambiguous_base_index(bp);
            const int j2 = unambiguous_base_index(bq);
            const int k2 = unambiguous_base_index(br);
            if (i2 >= 0 && j2 >= 0 && k2 >= 0) {
              pats.site_pattern[i2 * 16 + j2 * 4 + k2] += 1.0;
              continue;
            }

            const NucleotideDistribution dp = nucleotide_distribution(bp);
            const NucleotideDistribution dq = nucleotide_distribution(bq);
            const NucleotideDistribution dr = nucleotide_distribution(br);

            for (int i = 0; i < 4; ++i)
              for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                  pats.site_pattern[i * 16 + j * 4 + k] +=
                      dp[i] * dq[j] * dr[k];
          }
        }

        visit(std::move(pats));
      }
    }
  }
}

/**
 * Call visit(pats) for every k-tuple of taxa in a gene alignment, producing
 * one TuplePatterns at a time without intermediate storage.  k must be >= 1.
 *
 * Taxa within each TuplePatterns are in lexicographic order.
 * site_pattern has 4^k entries indexed as described in TuplePatterns.
 * When k=1, site_pattern[i] is the ACGT base count for taxon i.
 */
template <typename Visitor>
void for_each_ktuple_pattern(const GeneAlignment& gene,
                             size_t k,
                             Visitor&& visit) {
  if (k < 1)
    throw std::invalid_argument("for_each_ktuple_pattern: k must be >= 1");

  const size_t n = gene.num_taxa();
  const size_t len = gene.length();
  if (n < k)
    return;

  // Build a sorted index so taxa come out in lexicographic order.
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return gene.sequences[a].taxon < gene.sequences[b].taxon;
  });

  // Iterate over all (n choose k) combinations using an index vector.
  std::vector<size_t> combo(k);
  for (size_t i = 0; i < k; ++i)
    combo[i] = i;

  const bool use_compressed = !gene.compressed.empty();

  while (true) {
    std::vector<std::string> taxa(k);
    std::vector<size_t> seq_indices(k);
    for (size_t t = 0; t < k; ++t) {
      seq_indices[t] = order[combo[t]];
      taxa[t] = gene.sequences[seq_indices[t]].taxon;
    }

    TuplePatterns pats(std::move(taxa));
    pats.n_total = static_cast<double>(gene.length());

    std::vector<NucleotideDistribution> dists(k);

    if (use_compressed) {
      const auto& cs = gene.compressed;
      for (size_t u = 0; u < cs.n_patterns; ++u) {
        bool any_gap = false;
        bool all_unambiguous = true;
        size_t flat_idx = 0;
        if (!cs.all_defined[u]) {
          for (size_t t = 0; t < k; ++t) {
            if (is_gap(cs.base(u, seq_indices[t]))) {
              any_gap = true;
              break;
            }
          }
          if (any_gap)
            continue;
        }

        for (size_t t = 0; t < k; ++t) {
          const uint8_t b = cs.base(u, seq_indices[t]);
          const int bi = unambiguous_base_index(b);
          if (bi >= 0) {
            flat_idx = flat_idx * 4 + static_cast<size_t>(bi);
          } else {
            all_unambiguous = false;
            dists[t] = nucleotide_distribution(b);
          }
        }

        const double w = cs.weights[u];
        pats.n_counted += w;

        if (all_unambiguous) {
          pats.site_pattern[flat_idx] += w;
        } else {
          for (size_t t = 0; t < k; ++t) {
            const int bi = unambiguous_base_index(cs.base(u, seq_indices[t]));
            if (bi >= 0)
              dists[t] = nucleotide_distribution(cs.base(u, seq_indices[t]));
          }
          detail::accumulate_ktuple(pats.site_pattern, dists, 0, 0, w);
        }
      }
    } else {
      std::vector<const std::vector<uint8_t>*> seqs(k);
      for (size_t t = 0; t < k; ++t)
        seqs[t] = &gene.sequences[seq_indices[t]].bases;

      for (size_t s = 0; s < len; ++s) {
        bool any_gap = false;
        bool all_unambiguous = true;
        size_t flat_idx = 0;
        for (size_t t = 0; t < k; ++t) {
          const uint8_t b = (*seqs[t])[s];
          if (is_gap(b)) {
            any_gap = true;
            break;
          }
          const int bi = unambiguous_base_index(b);
          if (bi >= 0) {
            flat_idx = flat_idx * 4 + static_cast<size_t>(bi);
          } else {
            all_unambiguous = false;
            dists[t] = nucleotide_distribution(b);
          }
        }
        if (any_gap)
          continue;
        pats.n_counted += 1.0;

        if (all_unambiguous) {
          pats.site_pattern[flat_idx] += 1.0;
        } else {
          for (size_t t = 0; t < k; ++t) {
            const int bi = unambiguous_base_index((*seqs[t])[s]);
            if (bi >= 0)
              dists[t] = nucleotide_distribution((*seqs[t])[s]);
          }
          detail::accumulate_ktuple(pats.site_pattern, dists, 0, 0, 1.0);
        }
      }
    }

    visit(std::move(pats));

    // Advance to next combination in lexicographic order.
    int pos = static_cast<int>(k) - 1;
    while (pos >= 0 && combo[pos] == n - k + static_cast<size_t>(pos))
      --pos;
    if (pos < 0)
      break;
    ++combo[pos];
    for (size_t i = static_cast<size_t>(pos) + 1; i < k; ++i)
      combo[i] = combo[i - 1] + 1;
  }
}

/**
 * Compute site-pattern counts for a specific k-tuple of taxa in one gene.
 *
 * taxa may be in any order; they are sorted lexicographically internally.
 * k must be >= 1.  When k=1, site_pattern[i] is the ACGT base count
 * (0=A..3=T).  Throws std::invalid_argument if any taxon is missing or
 * names are duplicated.
 */
inline TuplePatterns site_patterns_for_ktaxa(const GeneAlignment& gene,
                                             std::vector<std::string> taxa) {
  std::sort(taxa.begin(), taxa.end());
  if (taxa.size() < 1)
    throw std::invalid_argument(
        "site_patterns_for_ktaxa: need at least 1 taxon");
  for (size_t t = 0; t + 1 < taxa.size(); ++t)
    if (taxa[t] == taxa[t + 1])
      throw std::invalid_argument(
          "site_patterns_for_ktaxa: duplicate taxon names");

  const size_t k = taxa.size();

  // Find the sequence index for each taxon.
  std::vector<size_t> seq_idx(k, SIZE_MAX);
  for (size_t i = 0; i < gene.sequences.size(); ++i)
    for (size_t t = 0; t < k; ++t)
      if (gene.sequences[i].taxon == taxa[t])
        seq_idx[t] = i;

  for (size_t t = 0; t < k; ++t)
    if (seq_idx[t] == SIZE_MAX)
      throw std::invalid_argument("site_patterns_for_ktaxa: taxon '" + taxa[t] +
                                  "' not found in alignment '" +
                                  gene.gene_name + "'");

  TuplePatterns pats(taxa);
  pats.n_total = static_cast<double>(gene.length());

  std::vector<NucleotideDistribution> dists(k);

  if (!gene.compressed.empty()) {
    const auto& cs = gene.compressed;
    for (size_t u = 0; u < cs.n_patterns; ++u) {
      bool any_gap = false;
      bool all_unambiguous = true;
      size_t flat_idx = 0;

      if (!cs.all_defined[u]) {
        for (size_t t = 0; t < k; ++t) {
          if (is_gap(cs.base(u, seq_idx[t]))) {
            any_gap = true;
            break;
          }
        }
        if (any_gap)
          continue;
      }

      for (size_t t = 0; t < k; ++t) {
        const uint8_t b = cs.base(u, seq_idx[t]);
        const int bi = unambiguous_base_index(b);
        if (bi >= 0) {
          flat_idx = flat_idx * 4 + static_cast<size_t>(bi);
        } else {
          all_unambiguous = false;
          dists[t] = nucleotide_distribution(b);
        }
      }

      const double w = cs.weights[u];
      pats.n_counted += w;

      if (all_unambiguous) {
        pats.site_pattern[flat_idx] += w;
      } else {
        for (size_t t = 0; t < k; ++t) {
          const int bi = unambiguous_base_index(cs.base(u, seq_idx[t]));
          if (bi >= 0)
            dists[t] = nucleotide_distribution(cs.base(u, seq_idx[t]));
        }
        detail::accumulate_ktuple(pats.site_pattern, dists, 0, 0, w);
      }
    }
  } else {
    const size_t len = gene.length();
    for (size_t s = 0; s < len; ++s) {
      bool any_gap = false;
      bool all_unambiguous = true;
      size_t flat_idx = 0;
      for (size_t t = 0; t < k; ++t) {
        const uint8_t b = gene.sequences[seq_idx[t]].bases[s];
        if (is_gap(b)) {
          any_gap = true;
          break;
        }
        const int bi = unambiguous_base_index(b);
        if (bi >= 0) {
          flat_idx = flat_idx * 4 + static_cast<size_t>(bi);
        } else {
          all_unambiguous = false;
          dists[t] = nucleotide_distribution(b);
        }
      }
      if (any_gap)
        continue;
      pats.n_counted += 1.0;

      if (all_unambiguous) {
        pats.site_pattern[flat_idx] += 1.0;
      } else {
        for (size_t t = 0; t < k; ++t) {
          const uint8_t b = gene.sequences[seq_idx[t]].bases[s];
          const int bi = unambiguous_base_index(b);
          if (bi >= 0)
            dists[t] = nucleotide_distribution(b);
        }
        detail::accumulate_ktuple(pats.site_pattern, dists, 0, 0, 1.0);
      }
    }
  }

  return pats;
}

/**
 * Compute aggregate site-pattern counts for a specific k-tuple of taxa
 * across all genes.
 * k must be >= 1.  When k=1, site_pattern[i] is the ACGT base count.
 * Throws if names are duplicated or any taxon is absent
 * from every gene.
 */
inline TuplePatterns aggregate_site_patterns(
    const MultiGeneAlignment& alignment,
    std::vector<std::string> taxa) {
  std::sort(taxa.begin(), taxa.end());
  if (taxa.size() < 1)
    throw std::invalid_argument(
        "aggregate_site_patterns: need at least 1 taxon");
  for (size_t t = 0; t + 1 < taxa.size(); ++t)
    if (taxa[t] == taxa[t + 1])
      throw std::invalid_argument(
          "aggregate_site_patterns: duplicate taxon names");

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
      throw std::invalid_argument("aggregate_site_patterns: taxon '" + taxa[t] +
                                  "' not found in any gene");
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
 * Compute aggregate site pattern counts by streaming over all genes.
 *
 * For each triplet, contributions are included only from genes that contain
 * all three taxa.  n_total and n_counted are summed accordingly.
 * Only one gene's worth of per-triplet data is live at a time.
 */
inline std::vector<TuplePatterns> compute_aggregate_site_patterns(
    const MultiGeneAlignment& alignment) {
  using Key = std::tuple<std::string, std::string, std::string>;
  std::map<Key, TuplePatterns> acc;
  int n_skipped = 0;

  for (const auto& gene : alignment.genes) {
    if (gene.num_taxa() < 3) {
      ++n_skipped;
      continue;
    }
    for_each_site_pattern(gene, [&](TuplePatterns pats) {
      Key key(pats.taxa[0], pats.taxa[1], pats.taxa[2]);
      TuplePatterns& a = acc[key];
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

  std::vector<TuplePatterns> result;
  result.reserve(acc.size());
  for (auto& kv : acc)
    result.push_back(std::move(kv.second));
  return result;
}

// ---------------------------------------------------------------------------
// Named public API
// ---------------------------------------------------------------------------

/**
 * Compute aggregate site-pattern counts for all n choose 3 taxon triplets
 * across all genes.  Equivalent to compute_aggregate_site_patterns().
 */
inline std::vector<TuplePatterns> triplet_site_patterns(
    const MultiGeneAlignment& alignment) {
  return compute_aggregate_site_patterns(alignment);
}

/**
 * Compute aggregate site-pattern counts for a specific taxon triplet across
 * all genes.
 */
inline TuplePatterns triplet_site_pattern(const MultiGeneAlignment& alignment,
                                          const std::string& taxon_a,
                                          const std::string& taxon_b,
                                          const std::string& taxon_c) {
  return aggregate_site_patterns(alignment, {taxon_a, taxon_b, taxon_c});
}

/**
 * Compute aggregate site-pattern counts for all n choose 4 taxon quartets
 * across all genes.  Only genes containing all four taxa contribute.
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
  for (auto& kv : acc)
    result.push_back(std::move(kv.second));
  return result;
}

/**
 * Compute aggregate site-pattern counts for a specific taxon quartet across
 * all genes.  Equivalent to aggregate_site_patterns() with k=4.
 */
inline TuplePatterns quartet_site_pattern(const MultiGeneAlignment& alignment,
                                          const std::string& taxon_a,
                                          const std::string& taxon_b,
                                          const std::string& taxon_c,
                                          const std::string& taxon_d) {
  return aggregate_site_patterns(alignment,
                                 {taxon_a, taxon_b, taxon_c, taxon_d});
}

/**
 * Compute aggregate site-pattern counts for a user-specified k-tuple of taxa
 * (k >= 1) across all genes.  When k=1, site_pattern[i] is the ACGT base
 * count (0=A..3=T).
 */
inline TuplePatterns site_pattern(const MultiGeneAlignment& alignment,
                                  std::vector<std::string> taxa) {
  return aggregate_site_patterns(alignment, std::move(taxa));
}

// ---------------------------------------------------------------------------
// Single-gene overloads (no aggregation, no skipping)
// ---------------------------------------------------------------------------

/** All n choose 3 triplet site-pattern counts for one gene. */
inline std::vector<TuplePatterns> triplet_site_patterns(
    const GeneAlignment& gene) {
  std::vector<TuplePatterns> result;
  for_each_site_pattern(
      gene, [&](TuplePatterns p) { result.push_back(std::move(p)); });
  return result;
}

/** Site-pattern counts for a specific triplet within one gene.
 *  Throws std::invalid_argument if any taxon is absent. */
inline TuplePatterns triplet_site_pattern(const GeneAlignment& gene,
                                          const std::string& taxon_a,
                                          const std::string& taxon_b,
                                          const std::string& taxon_c) {
  return site_patterns_for_ktaxa(gene, {taxon_a, taxon_b, taxon_c});
}

/** All n choose 4 quartet site-pattern counts for one gene. */
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

/** Site-pattern counts for a specific k-tuple (k >= 1) within one gene.
 *  When k=1, site_pattern[i] is the ACGT base count (0=A..3=T).
 *  Throws std::invalid_argument if any taxon is absent. */
inline TuplePatterns site_pattern(const GeneAlignment& gene,
                                  std::vector<std::string> taxa) {
  return site_patterns_for_ktaxa(gene, std::move(taxa));
}

// -- JC69 triplet weight convenience overloads --------------------------------
// (These depend on site_patterns_for_ktaxa / aggregate_site_patterns above,
//  so they live here rather than in triplet_weights.hpp.)

// Permute weights from sorted-taxa order back to caller-supplied order.
// sorted_weights[i] has outgroup sorted_taxa[2-i].
inline std::array<double, 3> permute_triplet_weights(
    const std::array<double, 3>& sorted_weights,
    const std::vector<std::string>& sorted_taxa,
    const std::string& taxon_a,
    const std::string& taxon_b,
    const std::string& taxon_c) {
  // Find each caller taxon's position in sorted order.
  // sorted_taxa has exactly 3 elements; each taxon matches exactly one.
  size_t pos_a = (taxon_a == sorted_taxa[0])   ? 0
                 : (taxon_a == sorted_taxa[1]) ? 1
                                               : 2;
  size_t pos_b = (taxon_b == sorted_taxa[0])   ? 0
                 : (taxon_b == sorted_taxa[1]) ? 1
                                               : 2;
  size_t pos_c = (taxon_c == sorted_taxa[0])   ? 0
                 : (taxon_c == sorted_taxa[1]) ? 1
                                               : 2;
  // Outgroup at sorted position j -> sorted_weights[2-j]
  return {sorted_weights[2 - pos_c],
          sorted_weights[2 - pos_b],
          sorted_weights[2 - pos_a]};
}

/**
 * Compute JC69 triplet weights from a multi-gene alignment.
 *
 * Results are in caller-supplied order:
 *   result[0] -- (taxon_a, taxon_b) | taxon_c
 *   result[1] -- (taxon_a, taxon_c) | taxon_b
 *   result[2] -- (taxon_b, taxon_c) | taxon_a
 */
inline std::array<double, 3> compute_jc69_triplet_weights(
    const MultiGeneAlignment& alignment,
    const std::string& taxon_a,
    const std::string& taxon_b,
    const std::string& taxon_c,
    double theta) {
  auto site_patterns =
      aggregate_site_patterns(alignment, {taxon_a, taxon_b, taxon_c});
  auto sorted_log_weights = compute_jc69_triplet_weights(site_patterns, theta);
  // const double scale = log_sum_exp(sorted_log_weights);
  // for (double& w : sorted_log_weights)
  //   w -= scale;
  std::array<double, 3> sorted_weights = {std::exp(sorted_log_weights[0]),
                                          std::exp(sorted_log_weights[1]),
                                          std::exp(sorted_log_weights[2])};
  return permute_triplet_weights(
      sorted_weights, site_patterns.taxa, taxon_a, taxon_b, taxon_c);
}

/**
 * Compute JC69 triplet weights from a single gene.
 *
 * Results are in caller-supplied order:
 *   result[0] -- (taxon_a, taxon_b) | taxon_c
 *   result[1] -- (taxon_a, taxon_c) | taxon_b
 *   result[2] -- (taxon_b, taxon_c) | taxon_a
 */
inline std::array<double, 3> compute_jc69_triplet_weights(
    const GeneAlignment& gene,
    const std::string& taxon_a,
    const std::string& taxon_b,
    const std::string& taxon_c,
    double theta) {
  auto site_patterns =
      site_patterns_for_ktaxa(gene, {taxon_a, taxon_b, taxon_c});
  auto sorted_log_weights = compute_jc69_triplet_weights(site_patterns, theta);
  const double scale = log_sum_exp(sorted_log_weights);
  if (scale == -std::numeric_limits<double>::infinity()) {
    // All triplets are rejected.
    return {0.0, 0.0, 0.0};
  }
  // for (double& w : sorted_log_weights)
  //   w -= scale;
  std::array<double, 3> sorted_weights = {std::exp(sorted_log_weights[0]),
                                          std::exp(sorted_log_weights[1]),
                                          std::exp(sorted_log_weights[2])};
  return permute_triplet_weights(
      sorted_weights, site_patterns.taxa, taxon_a, taxon_b, taxon_c);
}

/**
 * STUB: Compute weights for the three possible rooted triplet topologies for
 * a specific taxon triplet, based on site-pattern counts across all genes.
 *
 * Not yet implemented -- always throws std::runtime_error.
 */
inline std::array<double, 3> compute_triplet_weights(
    const MultiGeneAlignment& alignment,
    const std::string& taxon_a,
    const std::string& taxon_b,
    const std::string& taxon_c) {
  auto pis = alignment.nucleotide_frequencies();
  auto ss = std::array<double, 6>{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};  // JC69 model

  return {0.0, 0.0, 0.0};
}

/**
 * STUB: Compute weights for the three possible rooted triplet topologies for
 * a specific taxon triplet, based on site-pattern counts in a single gene.
 *
 * Not yet implemented -- always throws std::runtime_error.
 */
inline std::array<double, 3> compute_triplet_weights(
    const GeneAlignment& gene,
    const std::string& taxon_a,
    const std::string& taxon_b,
    const std::string& taxon_c) {
  throw std::runtime_error("compute_triplet_weights: not yet implemented");
}

}  // namespace phylotriplet