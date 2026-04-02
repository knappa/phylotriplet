#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace phylotriplet {

/**
 * Partition taxa using the spectral Laplacian method.
 *
 * When the adjacency graph is connected, returns two components produced by
 * a Fiedler-vector (sign) bipartition.  When the graph is disconnected,
 * returns one entry per connected component (there may be more than two).
 *
 * @return  Vector of components; each component is a vector of taxa names.
 *          Always contains at least two entries (the caller should treat a
 *          single-entry result as an error).
 */
inline std::vector<std::vector<std::string>> spectral_laplacian_partition(
    const Eigen::MatrixXd& adj_matrix,
    const std::vector<std::string>& taxa) {
  // basic checks
  const auto n = adj_matrix.rows();
  if (adj_matrix.cols() != n) {
    throw std::invalid_argument("adj_matrix must be square");
  }
  if (taxa.size() != static_cast<size_t>(n)) {
    throw std::invalid_argument("taxa size must match adjacency dimension");
  }

  // degree and Laplacian
  Eigen::VectorXd degree = adj_matrix.rowwise().sum();
  Eigen::MatrixXd L = Eigen::MatrixXd(degree.asDiagonal()) - adj_matrix;

  // eigen-decomposition (symmetric)
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(L);
  if (es.info() != Eigen::Success) {
    throw std::runtime_error("eigen decomposition failed");
  }
  Eigen::VectorXd evals = es.eigenvalues();
  Eigen::MatrixXd evecs = es.eigenvectors();

  // tolerances for "close to" checks
  const double tol = 1e-8 * std::max(1.0, evals.cwiseAbs().maxCoeff());

  // count zero (or near-zero) eigenvalues
  std::vector<int> zero_idx;
  for (int i = 0; i < evals.size(); ++i) {
    if (std::abs(evals[i]) <= tol)
      zero_idx.push_back(i);
  }

  if (zero_idx.size() > 1) {
    // Graph is disconnected: BFS to discover each connected component.
    std::vector<int> visited(n, 0);
    std::vector<std::vector<std::string>> components;
    for (int i = 0; i < n; ++i) {
      if (visited[i])
        continue;
      std::vector<std::string> comp;
      std::queue<int> q;
      q.push(i);
      visited[i] = 1;
      while (!q.empty()) {
        int u = q.front();
        q.pop();
        comp.push_back(taxa[u]);
        for (int v = 0; v < n; ++v) {
          if (!visited[v] && adj_matrix(u, v) != 0.0) {
            visited[v] = 1;
            q.push(v);
          }
        }
      }
      components.push_back(std::move(comp));
    }
    return components;
  }

  // single zero eigenvalue: use Fiedler vector (second smallest eigenvalue)
  int second_idx = (evals.size() >= 2) ? 1 : 0;
  double second_val = evals[second_idx];

  // find multiplicity around second smallest
  std::vector<int> special_idx;
  for (int i = 0; i < evals.size(); ++i) {
    if (std::abs(evals[i] - second_val) <= tol)
      special_idx.push_back(i);
  }

  Eigen::VectorXd component_vec(n);
  if (special_idx.size() > 1) {
    // build subspace matrix M = transpose(evecs[:, special_idx]) -> size k x n
    const int k = static_cast<int>(special_idx.size());
    Eigen::MatrixXd M(k, n);
    for (int r = 0; r < k; ++r) {
      M.row(r) = evecs.col(special_idx[r]).transpose();
    }
    // use SVD and take the right-singular-vector associated with smallest
    // singular value
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullV);
    Eigen::MatrixXd V = svd.matrixV();  // V is n x n
    component_vec = V.col(V.cols() - 1);
  } else {
    component_vec = evecs.col(second_idx);
  }

  // partition by sign, zeros go to smaller side, ties to negative side
  int pos_count = 0;
  int neg_count = 0;
  for (int i = 0; i < n; ++i) {
    if (component_vec[i] > 0)
      ++pos_count;
    else if (component_vec[i] < 0)
      ++neg_count;
  }
  bool prefer_nonneg_as_a = (pos_count < neg_count);
  std::vector<std::string> component_a;
  std::vector<std::string> component_b;
  for (int i = 0; i < n; ++i) {
    if (prefer_nonneg_as_a) {
      if (component_vec[i] >= 0)
        component_a.push_back(taxa[i]);
      else
        component_b.push_back(taxa[i]);
    } else {
      if (component_vec[i] > 0)
        component_a.push_back(taxa[i]);
      else
        component_b.push_back(taxa[i]);
    }
  }
  return {component_a, component_b};
}

/**
 * Same as spectral_laplacian_partition but treats the adjacency as unweighted:
 * any non-zero entry becomes 1.0.
 */
inline std::vector<std::vector<std::string>>
unweighted_spectral_laplacian_partition(const Eigen::MatrixXd& adj_matrix,
                                        const std::vector<std::string>& taxa) {
  Eigen::MatrixXd bin = (adj_matrix.array() != 0.0).cast<double>();
  return spectral_laplacian_partition(bin, taxa);
}

}  // namespace phylotriplet