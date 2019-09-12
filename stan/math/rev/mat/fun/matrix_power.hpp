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
  std::vector<Eigen::MatrixXd> Mds;

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
    Eigen::Map<Eigen::MatrixXd>(M_, rows_, cols_) = M.val();
    Eigen::Map<matrix_vi>(adjMRef_, rows_, cols_) = M.vi();
    Eigen::Map<matrix_d> Md(M_, rows_, cols_);
    Eigen::Map<matrix_d> resultd(result_, rows_, cols_);
    resultd = Eigen::MatrixXd::Identity(rows_, cols_);
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
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        Eigen::MatrixXd X = Eigen::MatrixXd::Zero(rows_, cols_);
        X(i, j) = 1.0;
        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(rows_, cols_);
        for (int nn = 0; nn < n; nn++) {
          S += Mds[nn] * X * Mds[n - nn - 1];
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
    * @param M A square matrix.
    * @param n Exponent.
    * @return nth power of M. M^n = M * ... * M.
    * @throw std::invalid_argument if the exponent is negative or the matrix is not
    * square.
    */
inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> matrix_power(
    const Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>& M, const int n) {
  check_square("matrix_power", "M", M);
  if (n < 0)
    invalid_argument("matrix_power", "n", n, "is ", ", but must be >= 0!");

  internal::matrix_product_vari<Eigen::Dynamic, Eigen::Dynamic>* baseVari
      = new internal::matrix_product_vari<Eigen::Dynamic, Eigen::Dynamic>(M, n);

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> result(M.rows(), M.cols());
  result.vi()
      = Eigen::Map<matrix_vi>(baseVari->adjResultRef_, M.rows(), M.cols());

  return result;
}

}  // namespace math
}  // namespace stan
#endif
