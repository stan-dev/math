#ifndef STAN_MATH_REV_MAT_FUN_DETERMINANT_HPP
#define STAN_MATH_REV_MAT_FUN_DETERMINANT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

namespace internal {
template <int R, int C>
class determinant_vari : public vari {
  int rows_;
  int cols_;
  double* A_;
  vari** adjARef_;

 public:
  explicit determinant_vari(const Eigen::Matrix<var, R, C>& A)
      : vari(determinant_vari_calc(A)),
        rows_(A.rows()),
        cols_(A.cols()),
        A_(reinterpret_cast<double*>(ChainableStack::instance_->memalloc_.alloc(
            sizeof(double) * A.rows() * A.cols()))),
        adjARef_(
            reinterpret_cast<vari**>(ChainableStack::instance_->memalloc_.alloc(
                sizeof(vari*) * A.rows() * A.cols()))) {
    Eigen::Map<Eigen::MatrixXd>(A_, rows_, cols_) = A.val();
    Eigen::Map<matrix_vi>(adjARef_, rows_, cols_) = A.vi();
  }
  static double determinant_vari_calc(const Eigen::Matrix<var, R, C>& A) {
    return A.val().determinant();
  }
  virtual void chain() {
    Eigen::Map<matrix_vi>(adjARef_, rows_, cols_).adj()
        += (adj_ * val_)
           * Eigen::Map<Eigen::MatrixXd>(A_, rows_, cols_)
                 .inverse()
                 .transpose();
  }
};
}  // namespace internal

template <int R, int C>
inline var determinant(const Eigen::Matrix<var, R, C>& m) {
  check_square("determinant", "m", m);
  return var(new internal::determinant_vari<R, C>(m));
}

}  // namespace math
}  // namespace stan
#endif
