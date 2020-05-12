#ifndef STAN_MATH_FWD_FUN_READ_FVAR_HPP
#define STAN_MATH_FWD_FUN_READ_FVAR_HPP

#include <stan/math/fwd/meta.hpp>

namespace stan {
namespace math {

/**
 * Functor for extracting the values and tangents from a matrix of fvar.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template <typename EigFvar, typename EigOut>
class read_fvar_functor {
  const EigFvar& var_mat;
  EigOut& val_mat;

 public:
  read_fvar_functor(const EigFvar& arg1, EigOut& arg2)
      : var_mat(arg1), val_mat(arg2) {}

  inline decltype(auto) operator()(Eigen::Index row, Eigen::Index col) const {
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).val_;
    return var_mat.coeffRef(row, col).d_;
  }

  inline decltype(auto) operator()(Eigen::Index index) const {
    val_mat.coeffRef(index) = var_mat.coeffRef(index).val_;
    return var_mat.coeffRef(index).d_;
  }
};

/**
 * Function applying the read_fvar_functor to extract the values
 * and tangets of a given fvar matrix into separate matrices.
 *
 * @tparam EigFvar type of the Eigen container of fvar.
 * @tparam EigOut type of the Eigen containers to copy to
 * @param[in] FvarMat Input Eigen container of fvar.
 * @param[in] ValMat Output Eigen container of values.
 * @param[in] DMat Output Eigen container of tangents.
 */
template <typename EigFvar, typename EigOut>
inline void read_fvar(const EigFvar& FvarMat, EigOut& ValMat, EigOut& DMat) {
  DMat = EigOut::NullaryExpr(
      FvarMat.rows(), FvarMat.cols(),
      read_fvar_functor<const EigFvar, EigOut>(FvarMat, ValMat));
}

}  // namespace math
}  // namespace stan
#endif
