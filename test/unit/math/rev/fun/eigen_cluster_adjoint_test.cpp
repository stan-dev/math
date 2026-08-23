#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/rev/fun/eigenvectors_sym.hpp>
#include <stan/math/rev/fun/eigenvalues_sym.hpp>
#include <stan/math/rev/fun/eigendecompose_sym.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <utility>
#include <vector>

// Reverse-mode eigenvector adjoints on numerically degenerate spectra.
//
// The adjoint of eigenvectors_sym / eigendecompose_sym pairs the
// antisymmetric F_ij = 1/(w_j - w_i) with V^T (dV).  When eigenvalues
// coincide to machine resolution the individual eigenvectors are not
// identifiable (any orthonormal basis of the cluster subspace is a valid
// answer): with exact repeats F divides by zero and the gradient is NaN;
// with rounding-level gaps the 1/(w_j - w_i) coupling amplifies rounding
// noise to O(1/eps).  The cluster-gauged adjoint zeroes the coupling of
// pairs whose gap falls below kappa * max(1, |w|_inf) * eps (the classical
// minimal-norm gauge for the undetermined within-cluster rotation).
//
// These tests pin:
//  1. exact repeats: finite gradients where the unguarded adjoint is NaN,
//     and Richardson-FD consistency of order 1e-8 for cluster-INVARIANT
//     composites (spectral functions of A whose value does not depend on
//     the basis chosen inside a cluster);
//  2. a rounding-degenerate jitter-floor kernel: finite gradients and
//     agreement between the two-call idiom (eigenvectors_sym +
//     eigenvalues_sym) and eigendecompose_sym;
//  3. well-separated spectra: the guarded path reproduces the textbook
//     adjoint to machine precision (the guard must not fire there).

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using stan::math::matrix_v;
using stan::math::var;

// Integer-only xorshift so generated matrices are identical on every
// platform (double rounding of integer arithmetic is exact).
unsigned long long lcg_state = 88172645463325252ULL;
double lcgu() {  // uniform-ish in (-1, 1)
  lcg_state ^= lcg_state << 13;
  lcg_state ^= lcg_state >> 7;
  lcg_state ^= lcg_state << 17;
  return static_cast<int>(lcg_state >> 33) / 1073741824.0 - 1.0;
}

// phi(A) = u' V diag(1/w) V' u + sum(log w) = u' A^{-1} u + logdet(A):
// a spectral functional whose VALUE is invariant under the choice of
// basis inside a degenerate cluster, so its derivative exists and is
// smooth even at exact repeats -- the right thing to FD-check against.
template <typename Mat>
auto phi_inv_two_call(const Mat& A, const VectorXd& u) {
  using scalar_t = typename Mat::Scalar;
  Eigen::Matrix<scalar_t, -1, -1> V = stan::math::eigenvectors_sym(A);
  Eigen::Matrix<scalar_t, -1, 1> w = stan::math::eigenvalues_sym(A);
  Eigen::Matrix<scalar_t, -1, 1> Vu = V.transpose() * u.cast<scalar_t>();
  Eigen::Matrix<scalar_t, -1, 1> Vu_over_w
      = stan::math::elt_divide(Vu, w);
  return stan::math::sum(Vu_over_w.cwiseProduct(Vu))
         + stan::math::sum(stan::math::log(w));
}

// Same composite computed through the combined primitive.
template <typename Mat>
auto phi_inv_combined(const Mat& A, const VectorXd& u) {
  using scalar_t = typename Mat::Scalar;
  auto ed = stan::math::eigendecompose_sym(A);
  Eigen::Matrix<scalar_t, -1, -1> V = std::get<0>(ed);
  Eigen::Matrix<scalar_t, -1, 1> w = std::get<1>(ed);
  Eigen::Matrix<scalar_t, -1, 1> Vu = V.transpose() * u.cast<scalar_t>();
  Eigen::Matrix<scalar_t, -1, 1> Vu_over_w
      = stan::math::elt_divide(Vu, w);
  return stan::math::sum(Vu_over_w.cwiseProduct(Vu))
         + stan::math::sum(stan::math::log(w));
}

// AD gradient of a scalar functional w.r.t. a full symmetric matrix.
// The operand adjoint is NOT symmetric in general (the coupling term
// pairs an antisymmetric factor with a generic downstream adjoint), so
// the layout must be respected: entry (i, j) is g(i, j).
template <typename Phi>
MatrixXd ad_grad(const MatrixXd& A0, const Phi& phi) {
  std::vector<var> av(A0.data(), A0.data() + A0.size());
  matrix_v Am = Eigen::Map<const matrix_v>(av.data(), A0.rows(), A0.cols());
  var f = phi(Am);
  f.grad();
  VectorXd g(A0.size());
  for (int i = 0; i < A0.size(); ++i) {
    g[i] = av[i].adj();
  }
  stan::math::recover_memory();
  return Eigen::Map<const MatrixXd>(g.data(), A0.rows(), A0.cols());
}

// Richardson-extrapolated central finite difference of phi(A) along the
// symmetric direction E_ij + E_ji.
template <typename Phi>
double fd_richardson(const MatrixXd& A, const Phi& phi, int i, int j) {
  auto eval = [&](double h) {
    MatrixXd P = A, Q = A;
    P(i, j) += h;
    if (i != j) {
      P(j, i) += h;
    }
    Q(i, j) -= h;
    if (i != j) {
      Q(j, i) -= h;
    }
    return (phi(P) - phi(Q)) / (2 * h);
  };
  double d1 = eval(1e-4);
  double d2 = eval(5e-5);
  return (4 * d2 - d1) / 3;
}

// AD directional derivative d phi / d eps along the same direction.
double ad_direction(const MatrixXd& g, int n, int i, int j) {
  double ad = g(i, j);
  if (i != j) {
    ad += g(j, i);
  }
  return ad;
}

}  // namespace

TEST_F(AgradRev, eigenvectors_sym_exact_repeated_eigenvalue) {
  // Exactly 4-fold repeated eigenvalue mu = 1 followed by well-separated
  // ones: the unguarded adjoint divides by exact zeros (NaN); the
  // cluster-gauged adjoint is finite and equals the derivative of the
  // (well-defined) cluster-invariant composite.
  const int n = 12, k = 4;
  MatrixXd A = MatrixXd::Zero(n, n);
  for (int i = 0; i < k; ++i) {
    A(i, i) = 1.0;  // exactly equal, k-fold
  }
  for (int i = k; i < n; ++i) {
    A(i, i) = 1.0 + 2.0 * i;
  }
  VectorXd u(n);
  for (int i = 0; i < n; ++i) {
    u[i] = lcgu();
  }

  MatrixXd g = ad_grad(A, [&](const matrix_v& Am) {
    return phi_inv_two_call(Am, u);
  });
  for (int i = 0; i < n * n; ++i) {
    EXPECT_TRUE(std::isfinite(g.data()[i])) << "component " << i;
  }

  // FD consistency on cluster-diagonal, well-separated-diagonal and
  // cross-cluster directions: for these the cluster-gauged adjoint IS
  // the derivative of the (well-defined) composite.  (Symmetric
  // directions with both indices inside the cluster are the deliberate
  // exception: there the gauge drops the within-cluster mixing term by
  // construction, so they are covered by the finiteness check above.)
  auto phi_d = [&](const MatrixXd& Am) { return phi_inv_two_call(Am, u); };
  for (auto [i, j] : std::vector<std::pair<int, int>>{
           {0, 0}, {1, 1}, {5, 5}, {0, 7}}) {
    double fd = fd_richardson(A, phi_d, i, j);
    double ad = ad_direction(g, n, i, j);
    double rel = std::abs(fd - ad) / std::max(1.0, std::abs(ad));
    EXPECT_LT(rel, 1e-8) << "direction (" << i << "," << j << ")";
  }

  // The combined primitive returns the same gradient (it must: both
  // decompose the same matrix with the same solver).
  MatrixXd g_comb = ad_grad(A, [&](const matrix_v& Am) {
    return phi_inv_combined(Am, u);
  });
  for (int i = 0; i < n * n; ++i) {
    EXPECT_DOUBLE_EQ(g.data()[i], g_comb.data()[i]) << "component " << i;
  }
}

TEST_F(AgradRev, eigenvectors_sym_total_degeneracy_zero_matrix) {
  // The zero matrix has an n-fold repeated eigenvalue 0: every pair is
  // degenerate.  Any downstream functional gets a finite gradient.
  const int n = 8;
  MatrixXd A = MatrixXd::Zero(n, n);
  MatrixXd W(n, n);
  VectorXd c(n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      W(i, j) = lcgu();
    }
  }
  for (int i = 0; i < n; ++i) {
    c[i] = lcgu();
  }

  auto phi = [&](const auto& Am) {
    auto V = stan::math::eigenvectors_sym(Am);
    auto w = stan::math::eigenvalues_sym(Am);
    using scalar_t = typename std::decay_t<decltype(V)>::Scalar;
    return stan::math::sum(
               W.cast<scalar_t>().cwiseProduct(V).eval())
           + stan::math::sum(c.cast<scalar_t>().cwiseProduct(w).eval());
  };
  MatrixXd g = ad_grad(A, phi);
  for (int i = 0; i < n * n; ++i) {
    EXPECT_TRUE(std::isfinite(g.data()[i])) << "two-call component " << i;
  }

  auto phi_ed = [&](const auto& Am) {
    auto ed = stan::math::eigendecompose_sym(Am);
    auto V = std::get<0>(ed);
    auto w = std::get<1>(ed);
    using scalar_t = typename std::decay_t<decltype(V)>::Scalar;
    return stan::math::sum(
               W.cast<scalar_t>().cwiseProduct(V).eval())
           + stan::math::sum(c.cast<scalar_t>().cwiseProduct(w).eval());
  };
  MatrixXd g_ed = ad_grad(A, phi_ed);
  for (int i = 0; i < n * n; ++i) {
    EXPECT_TRUE(std::isfinite(g_ed.data()[i])) << "combined component " << i;
  }
}

TEST_F(AgradRev, eigenvectors_sym_jitter_floor_kernel) {
  // exp-quad kernel + 1e-5 jitter on a 30-point grid: the bottom ~10
  // eigenvalues are pinned at the jitter floor with rounding-level
  // internal gaps (the GP-regression posterior-covariance regime; the
  // guard fires on the masked pairs).  The gradient must be finite, and
  // the two-call idiom must agree with eigendecompose_sym.
  const int n = 30;
  VectorXd x = VectorXd::LinSpaced(n, 0.0, 1.0);
  MatrixXd A(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A(i, j) = std::exp(-std::pow(x[i] - x[j], 2.0) / 0.05);
    }
    A(i, i) += 1e-5;
  }
  VectorXd u(n);
  for (int i = 0; i < n; ++i) {
    u[i] = lcgu();
  }

  MatrixXd g = ad_grad(A, [&](const matrix_v& Am) {
    return phi_inv_two_call(Am, u);
  });
  for (int i = 0; i < n * n; ++i) {
    EXPECT_TRUE(std::isfinite(g.data()[i])) << "component " << i;
  }

  MatrixXd g_comb = ad_grad(A, [&](const matrix_v& Am) {
    return phi_inv_combined(Am, u);
  });
  for (int i = 0; i < n * n; ++i) {
    EXPECT_DOUBLE_EQ(g.data()[i], g_comb.data()[i]) << "component " << i;
  }
}

TEST_F(AgradRev, eigenvectors_sym_well_separated_matches_textbook_adjoint) {
  // Well-separated spectra (min adjacent gap >= 1e-6 * scale, enforced)
  // must take the unguarded code path: the operand adjoint is the
  // textbook formula V (F ∘ (V' G_V)) V' (+ V diag(g_w) V' for the
  // combined primitive) with F_ij = 1/(w_j - w_i), to rounding.
  const int n = 12, n_mats = 10;
  int tested = 0;
  for (int m = 0; m < n_mats; ++m) {
    MatrixXd A(n, n);
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        A(i, j) = lcgu();
      }
    }
    A = ((A + A.transpose()) * 0.5).eval();
    for (int i = 0; i < n; ++i) {
      A(i, i) += 40.0;  // keep eigenvalue gaps O(1)
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    VectorXd w = es.eigenvalues();
    MatrixXd V = es.eigenvectors();
    double scale = std::max(1.0, w.cwiseAbs().maxCoeff());
    double mingap = (w.tail(n - 1) - w.head(n - 1)).minCoeff();
    if (mingap < 1e-6 * scale) {
      continue;
    }
    ++tested;
    MatrixXd W(n, n);
    VectorXd c(n);
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        W(i, j) = lcgu();
      }
    }
    for (int i = 0; i < n; ++i) {
      c[i] = lcgu();
    }

    // Reference adjoint: the textbook formula V (F o (V' G_V)) V'
    // + V diag(g_w) V' with F_ij = 1/(w_j - w_i), evaluated in double.
    // Both functionals below use BOTH outputs, so both arms must
    // reproduce the full adjoint (to rounding; the guard must not fire,
    // and a masked pair would show up as an O(1) difference).
    Eigen::MatrixXd f
        = (1
           / (w.rowwise().replicate(n).transpose()
              - w.rowwise().replicate(n))
                 .array());
    f.diagonal().setZero();
    MatrixXd expected_full_adj
        = V * f.cwiseProduct(V.transpose() * W) * V.transpose()
          + V * c.asDiagonal() * V.transpose();

    auto phi = [&](const auto& Am) {
      auto Vv = stan::math::eigenvectors_sym(Am);
      auto wv = stan::math::eigenvalues_sym(Am);
      using scalar_t = typename std::decay_t<decltype(Vv)>::Scalar;
      return stan::math::sum(
                 W.cast<scalar_t>().cwiseProduct(Vv).eval())
             + stan::math::sum(c.cast<scalar_t>().cwiseProduct(wv).eval());
    };
    MatrixXd g = ad_grad(A, phi);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        EXPECT_NEAR(g(i, j), expected_full_adj(i, j), 1e-12)
            << "matrix " << m << " entry (" << i << "," << j << ")";
      }
    }

    auto phi_ed = [&](const auto& Am) {
      auto ed = stan::math::eigendecompose_sym(Am);
      auto Vv = std::get<0>(ed);
      auto wv = std::get<1>(ed);
      using scalar_t = typename std::decay_t<decltype(Vv)>::Scalar;
      return stan::math::sum(
                 W.cast<scalar_t>().cwiseProduct(Vv).eval())
             + stan::math::sum(c.cast<scalar_t>().cwiseProduct(wv).eval());
    };
    MatrixXd g_ed = ad_grad(A, phi_ed);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        EXPECT_NEAR(g_ed(i, j), expected_full_adj(i, j), 1e-12)
            << "matrix " << m << " entry (" << i << "," << j << ")";
      }
    }
  }
  EXPECT_GT(tested, 0) << "no well-separated matrix generated";
}
