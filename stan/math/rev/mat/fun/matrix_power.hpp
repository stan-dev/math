#ifndef STAN_MATH_REV_MAT_FUN_MATRIX_POWER_HPP
#define STAN_MATH_REV_MAT_FUN_MATRIX_POWER_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <vector>
#include <iostream>

namespace stan {
namespace math {

namespace internal {

template <int R, int C>
class matrix_product_vari : public vari {
 public:
  int rows_;
  int cols_;
  double* M_;
  double* result_;
  vari** adjMRef_;
  vari** adjResultRef_;
  std::vector<Eigen::Matrix<double, R, C> > Mds;

  explicit matrix_product_vari(const Eigen::Matrix<var, R, C>& M, const int n)
      : vari(0.0),
        rows_(M.rows()),
        cols_(M.cols()),
        M_(reinterpret_cast<double*>(ChainableStack::instance_->memalloc_.alloc(
            sizeof(double) * M.rows() * M.cols()))),
        result_(reinterpret_cast<double*>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * M.rows()
                                                       * M.cols()))),
        adjMRef_(
            reinterpret_cast<vari**>(ChainableStack::instance_->memalloc_.alloc(
                sizeof(vari*) * M.rows() * M.cols()))),
        adjResultRef_(
            reinterpret_cast<vari**>(ChainableStack::instance_->memalloc_.alloc(
                sizeof(vari*) * M.rows() * M.cols()))) {
    Eigen::Map<Eigen::Matrix<double, R, C> >(M_, rows_, cols_) = M.val();
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = M.vi();
    Eigen::Map<matrix_d> Md(M_, rows_, cols_);
    Eigen::Map<matrix_d> resultd(result_, rows_, cols_);
    resultd = Eigen::Matrix<double, R, C>::Identity(rows_, cols_);
    for (int i = 0; i < n; i++) {
      Mds.push_back(resultd);
      resultd *= Md;
    }
    Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
        = resultd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    int n = Mds.size();
    if (n == 0)
      return;
    if (n == 1) {
      Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_).adj()
          = Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_).adj();
      return;
    }
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        double dijdkl = 0.0;
        Eigen::Matrix<double, R, C> S
            = Eigen::Matrix<double, R, C>::Zero(rows_, cols_);
        for (int l = 0; l < cols_; l++) {
          S(i, l) += Mds[n - 1](j, l);
        }
        for (int nn = 1; nn < n - 1; nn++) {
          for (int k = 0; k < rows_; k++) {
            for (int l = 0; l < cols_; l++) {
              S(k, l) += Mds[nn](k, i) * Mds[n - nn - 1](j, l);
            }
          }
        }
        for (int k = 0; k < rows_; k++) {
          S(k, j) += Mds[n - 1](k, i);
        }
        Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_).adj()(i, j)
            += (S.array()
                * Eigen::Map<matrix_vi>(adjResultRef_, rows_, cols_)
                      .adj()
                      .array())
                   .sum();
      }
    }
  }
};

}  // namespace internal

/**
 * Returns the nth power of the specific matrix.
 *
 * @tparam R Number of rows in matrix.
 * @tparam C Number of columns in matrix.
 * @param M A square matrix.
 * @param n Exponent.
 * @return nth power of M. M^n = M * ... * M.
 * @throw std::domain_error if the matrix contains NaNs or infinities.
 * @throw std::invalid_argument if the exponent is negative or the matrix is not
 * square.
 */
template <int R, int C>
inline Eigen::Matrix<var, R, C> matrix_power(const Eigen::Matrix<var, R, C>& M,
                                             const int n) {
  check_square("matrix_power", "M", M);
  if (n < 0)
    invalid_argument("matrix_power", "n", n, "is ", ", but must be >= 0!");
  if (M.rows() == 0)
    invalid_argument("matrix_power", "M.rows()", M.rows(), "is ",
                     ", but must be > 0!");
  check_finite("matrix_power", "M", M);
  internal::matrix_product_vari<R, C>* baseVari
      = new internal::matrix_product_vari<R, C>(M, n);
  Eigen::Matrix<var, R, C> result(M.rows(), M.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, M.rows(), M.cols());
  return result;
}

}  // namespace math
}  // namespace stan
#endif
