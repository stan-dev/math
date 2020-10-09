#ifndef STAN_MATH_REV_FUN_MATRIX_POWER_HPP
#define STAN_MATH_REV_FUN_MATRIX_POWER_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <int R, int C>
class matrix_product_vari_n0 : public vari {
 public:
  vari** adjMnRef_;

  explicit matrix_product_vari_n0(const Eigen::Matrix<var, R, C>& M)
      : vari(0.0),
        adjMnRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            M.rows() * M.cols())) {
    matrix_d Mnd = Eigen::Matrix<double, R, C>::Identity(M.rows(), M.cols());
    Eigen::Map<matrix_vi>(adjMnRef_, M.rows(), M.cols())
        = Mnd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() { return; }
};

template <int R, int C>
class matrix_product_vari_n1 : public vari {
 public:
  int rows_;
  int cols_;
  vari** adjMRef_;
  vari** adjMnRef_;

  explicit matrix_product_vari_n1(const Eigen::Matrix<var, R, C>& M)
      : vari(0.0),
        rows_(M.rows()),
        cols_(M.cols()),
        adjMRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            M.rows() * M.cols())),
        adjMnRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            M.rows() * M.cols())) {
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = M.vi();
    Eigen::Map<matrix_vi>(adjMnRef_, M.rows(), M.cols())
        = M.val().unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_).adj()
        = Eigen::Map<matrix_vi>(adjMnRef_, rows_, cols_).adj();
  }
};

template <int R, int C>
class matrix_product_vari : public vari, public chainable_alloc {
 public:
  int rows_;
  int cols_;
  vari** adjMRef_;
  vari** adjMnRef_;
  std::vector<Eigen::Matrix<double, R, C> > Mds;  // M, M^2, ..., M^(n-1).
  int n;

  explicit matrix_product_vari(const Eigen::Matrix<var, R, C>& M, const int n)
      : vari(0.0),
        rows_(M.rows()),
        cols_(M.cols()),
        adjMRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            M.rows() * M.cols())),
        adjMnRef_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            M.rows() * M.cols())),
        n(n) {
    Mds.resize(n - 1);
    Mds[0] = M.val();
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = M.vi();
    for (int i = 1; i < n - 1; i++) {
      Mds[i] = Mds[i - 1] * Mds[0];
    }
    Eigen::Matrix<double, R, C> Mnd = Mds[n - 2] * Mds[0];
    Eigen::Map<matrix_vi>(adjMnRef_, rows_, cols_)
        = Mnd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    for (int j = 0; j < cols_; j++) {
      for (int i = 0; i < rows_; i++) {
        double dijdkl = 0.0;
        Eigen::Matrix<double, R, C> S
            = Eigen::Matrix<double, R, C>::Zero(rows_, cols_);
        for (int l = 0; l < cols_; l++) {
          S(i, l) += Mds[n - 2](j, l);
        }
        for (int nn = 1; nn < n - 1; nn++) {
          for (int l = 0; l < cols_; l++) {
            for (int k = 0; k < rows_; k++) {
              S(k, l) += Mds[nn - 1](k, i) * Mds[n - nn - 2](j, l);
            }
          }
        }
        for (int k = 0; k < rows_; k++) {
          S(k, j) += Mds[n - 2](k, i);
        }
        Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_).adj()(i, j)
            += (S.array()
                * Eigen::Map<matrix_vi>(adjMnRef_, rows_, cols_).adj().array())
                   .sum();
      }
    }
  }
};

}  // namespace internal

/**
 * Returns the nth power of the specific matrix. M^n = M * M * ... * M.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param[in] M a square matrix
 * @param[in] n exponent
 * @return nth power of M
 * @throw std::domain_error if the matrix contains NaNs or infinities.
 * @throw std::invalid_argument if the exponent is negative or the matrix is not
 * square.
 */
template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
matrix_power(const EigMat& M, const int n) {
  using T = value_type_t<EigMat>;
  constexpr int R = EigMat::RowsAtCompileTime;
  constexpr int C = EigMat::ColsAtCompileTime;

  check_square("matrix_power", "M", M);
  if (n < 0)
    invalid_argument("matrix_power", "n", n, "is ", ", but must be >= 0!");
  if (M.rows() == 0)
    invalid_argument("matrix_power", "M.rows()", M.rows(), "is ",
                     ", but must be > 0!");
  const auto& M_ref = to_ref(M);
  check_finite("matrix_power", "M", M_ref);
  if (n == 0) {
    auto* baseVari = new internal::matrix_product_vari_n0<R, C>(M_ref);
    Eigen::Matrix<var, R, C> Mn(M.rows(), M.cols());
    Mn.vi() = Eigen::Map<matrix_vi>(baseVari->adjMnRef_, M.rows(), M.cols());
    return Mn;
  } else if (n == 1) {
    auto* baseVari = new internal::matrix_product_vari_n1<R, C>(M_ref);
    Eigen::Matrix<var, R, C> Mn(M.rows(), M.cols());
    Mn.vi() = Eigen::Map<matrix_vi>(baseVari->adjMnRef_, M.rows(), M.cols());
    return Mn;
  } else {
    auto* baseVari = new internal::matrix_product_vari<R, C>(M_ref, n);
    Eigen::Matrix<var, R, C> Mn(M.rows(), M.cols());
    Mn.vi() = Eigen::Map<matrix_vi>(baseVari->adjMnRef_, M.rows(), M.cols());
    return Mn;
  }
}

}  // namespace math
}  // namespace stan
#endif
