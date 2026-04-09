#define CATCH_CONFIG_MAIN
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

#include "catch.hpp"
#include "phylotriplet/optimize.hpp"

using namespace phylotriplet;

// Helper: build a TuplePatterns from known site pattern counts for 3 taxa.
static TuplePatterns make_patterns(const std::array<std::string, 3>& taxa,
                                   const std::array<double, 64>& counts) {
	std::vector<std::string> sorted_taxa(taxa.begin(), taxa.end());
	TuplePatterns tp(sorted_taxa);
	double total = 0.0;
	for (size_t i = 0; i < 64; ++i) {
		tp.site_pattern[i] = counts[i];
		total += counts[i];
	}
	tp.n_total = total;
	tp.n_counted = total;
	return tp;
}

// Generate expected site patterns from a known JC69 tree.
// Branch lengths in standard JC69 parameterization (rate-1 Q).
static std::array<double, 64> jc69_expected_patterns(double t_ai, double t_bi,
                                                     double t_cr,
                                                     double t_ir) {
	auto p = [](double t, int from, int to) -> double {
		double e = std::exp(-4.0 * t / 3.0);
		return (from == to) ? (0.25 + 0.75 * e) : (0.25 - 0.25 * e);
	};

	std::array<double, 64> patterns;
	patterns.fill(0.0);
	for (int a = 0; a < 4; ++a)
		for (int b = 0; b < 4; ++b)
			for (int c = 0; c < 4; ++c) {
				double prob = 0.0;
				for (int x = 0; x < 4; ++x)
					for (int y = 0; y < 4; ++y)
						prob += 0.25 * p(t_ir, y, x) * p(t_ai, x, a) *
						        p(t_bi, x, b) * p(t_cr, y, c);
				patterns[site_pattern_index(a, b, c)] = prob;
			}
	return patterns;
}

// --- JC69 constrained optimizer tests ----------------------------------------

TEST_CASE("extract_jc69_class_counts: proportions sum to 1") {
	auto patterns = jc69_expected_patterns(0.1, 0.2, 0.3, 0.15);
	const double n_sites = 10000.0;
	for (double& p : patterns) p *= n_sites;
	auto tp = make_patterns({"A", "B", "C"}, patterns);

	auto c = extract_jc69_class_counts(tp, 0, 1, 2);
	double total = c.n_xxx + c.n_xxy + c.n_xyx + c.n_xyz;
	REQUIRE(total == Approx(1.0).margin(1e-10));
}

TEST_CASE("optimize_jc69_topology_fixed_theta: recovers tau1 < tau2 for "
          "correct topology") {
	// Tree: ((A:0.05, B:0.05):0.2, C:0.3)
	// tau1 (cherry depth) < tau2 (root depth)
	auto patterns = jc69_expected_patterns(0.05, 0.05, 0.3, 0.2);
	const double n_sites = 100000.0;
	for (double& p : patterns) p *= n_sites;
	auto tp = make_patterns({"A", "B", "C"}, patterns);

	auto counts = extract_jc69_class_counts(tp, 0, 1, 2);  // (A,B)|C

	std::array<double, 2> params;
	double ll = optimize_jc69_topology_fixed_theta(counts, 0.01, &params);

	INFO("tau1=" << params[0] << " tau2=" << params[1]);
	REQUIRE(ll < 0.0);
	REQUIRE(params[0] >= 0.0);
	REQUIRE(params[1] >= params[0]);  // tau2 >= tau1 by construction (delta >= 0)
	REQUIRE(params[1] > params[0] + 0.001);  // tau2 strictly > tau1 for this tree
}

TEST_CASE("compute_jc69_triplet_weights_optimized: correct topology wins") {
	// Tree: ((A:0.05, B:0.05):0.3, C:0.5) -- strong cherry signal
	auto patterns = jc69_expected_patterns(0.05, 0.05, 0.5, 0.3);
	const double n_sites = 100000.0;
	for (double& p : patterns) p *= n_sites;
	auto tp = make_patterns({"A", "B", "C"}, patterns);

	auto weights = compute_jc69_triplet_weights_optimized(tp, 0.01);

	INFO("Weights: " << weights[0] << ", " << weights[1] << ", " << weights[2]);

	// Topology 0 = (A,B)|C should have highest log-likelihood
	REQUIRE(weights[0] > weights[1]);
	REQUIRE(weights[0] > weights[2]);
}

TEST_CASE("compute_jc69_triplet_weights_optimized: topology (B,C)|A wins") {
	// Tree: ((B:0.05, C:0.05):0.3, A:0.5)
	// In our notation: cherry=(B,C), outgroup=A.
	// Site patterns generated with: a=A outgroup, i connects B and C.
	// But jc69_expected_patterns takes (t_ai, t_bi, t_cr, t_ir) for tree
	//   ((taxa[0]:t_ai, taxa[1]:t_bi):t_ir, taxa[2]:t_cr)
	// We want ((B,C)|A), so we generate with taxa = {A, B, C} and
	// the tree ((A:pendant, B:0.05):0.3, C:0.05)... no, that's wrong.
	//
	// The easiest way: generate patterns for ((A:0.05, B:0.05):0.3, C:0.5)
	// but name the taxa {B, C, A} so that the cherry is (B,C).
	// Since site_patterns are indexed by sorted taxa, generate for {A,B,C}
	// where the TRUE tree has cherry=(B,C) = (taxa[1],taxa[2]).
	// This means: ((B:0.05, C:0.05):0.3, A:0.5)
	// Rewrite as: taxa[0]=A is outgroup, taxa[1]=B and taxa[2]=C are cherry.
	// In jc69_expected_patterns format (positions 0,1,2 = A,B,C):
	//   A is outgroup -> connected to root with t_cr (but A is position 0)
	// This doesn't map cleanly to jc69_expected_patterns which assumes
	// positions 0,1 are the cherry. Let's just swap the pattern indices.

	// Generate with cherry = (pos0, pos1), outgroup = pos2
	auto base = jc69_expected_patterns(0.05, 0.05, 0.5, 0.3);
	// Now reindex so cherry is (pos1, pos2) and outgroup is pos0.
	// Original: index = pos0*16 + pos1*4 + pos2
	// We want: the cherry pair to be (B=pos1, C=pos2), outgroup A=pos0.
	// Swap: new_pattern[a*16 + b*4 + c] = base[b*16 + c*4 + a]
	//   where in base, positions 0,1 are cherry and 2 is outgroup.
	std::array<double, 64> patterns;
	patterns.fill(0.0);
	for (int a = 0; a < 4; ++a)
		for (int b = 0; b < 4; ++b)
			for (int c = 0; c < 4; ++c)
				patterns[a * 16 + b * 4 + c] = base[b * 16 + c * 4 + a];

	const double n_sites = 100000.0;
	for (double& p : patterns) p *= n_sites;
	auto tp = make_patterns({"A", "B", "C"}, patterns);

	auto weights = compute_jc69_triplet_weights_optimized(tp, 0.01);

	INFO("Weights: " << weights[0] << ", " << weights[1] << ", " << weights[2]);

	// Topology 2 = (B,C)|A should win
	REQUIRE(weights[2] > weights[0]);
	REQUIRE(weights[2] > weights[1]);
}

TEST_CASE("compute_jc69_triplet_weights_optimized: star tree gives similar "
          "weights") {
	// Star tree: internal branch = 0
	auto patterns = jc69_expected_patterns(0.2, 0.2, 0.2, 0.0);
	const double n_sites = 100000.0;
	for (double& p : patterns) p *= n_sites;
	auto tp = make_patterns({"A", "B", "C"}, patterns);

	auto weights = compute_jc69_triplet_weights_optimized(tp, 0.01);

	INFO("Weights: " << weights[0] << ", " << weights[1] << ", " << weights[2]);

	// All topologies should yield similar log-likelihoods
	REQUIRE(weights[0] == Approx(weights[1]).margin(0.01));
	REQUIRE(weights[0] == Approx(weights[2]).margin(0.01));
}

// --- GTR optimizer tests (branch-length recovery) ----------------------------

TEST_CASE("optimize_triplet_likelihood: recovers pendant branch lengths") {
	const double t_ai = 0.1, t_bi = 0.2, t_cr = 0.3, t_ir = 0.15;
	const double scale = 4.0 / 3.0;

	auto patterns = jc69_expected_patterns(t_ai, t_bi, t_cr, t_ir);
	const double n_sites = 100000.0;
	for (double& p : patterns) p *= n_sites;
	auto tp = make_patterns({"A", "B", "C"}, patterns);

	const std::array<double, 4> pis = {{0.25, 0.25, 0.25, 0.25}};
	const std::array<double, 6> ss = {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
	const std::array<int, 3> taxon_order = {{0, 1, 2}};

	std::array<double, 4> best_bl;
	double ll = optimize_triplet_likelihood(pis, ss, tp, taxon_order, best_bl);

	INFO("Branch lengths: " << best_bl[0] << ", " << best_bl[1] << ", "
	                         << best_bl[2] << ", " << best_bl[3]);

	// Pendant branches are identifiable
	REQUIRE(best_bl[0] == Approx(t_ai * scale).margin(0.02));
	REQUIRE(best_bl[1] == Approx(t_bi * scale).margin(0.02));

	// Under JC69, root placement is unidentifiable; only total internal
	// path length (c->r + i->r) is constrained
	double total_internal = best_bl[2] + best_bl[3];
	double expected_internal = (t_cr + t_ir) * scale;
	REQUIRE(total_internal == Approx(expected_internal).margin(0.02));

	REQUIRE(ll < 0.0);
}
