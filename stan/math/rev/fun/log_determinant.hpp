#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline var log_determinant(const T& m) {
  check_square("determinant", "m", m);
  if (m.size() == 0) {
    return 0.0;
  }

  Eigen::FullPivHouseholderQR<promote_scalar_t<double, plain_type_t<T>>> hh
      = m.val().fullPivHouseholderQr();

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arena_m = m;
  arena_matrix<Eigen::MatrixXd> arena_hh_inv_t = hh.inverse().transpose();

  var log_det = hh.logAbsDeterminant();

  reverse_pass_callback([arena_m, log_det, arena_hh_inv_t]() mutable {
    arena_m.adj() += log_det.adj() * arena_hh_inv_t;
  });

  return log_det;
}

}  // namespace math
}  // namespace stan
#endif
