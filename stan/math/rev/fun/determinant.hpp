#ifndef STAN_MATH_REV_FUN_DETERMINANT_HPP
#define STAN_MATH_REV_FUN_DETERMINANT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline var determinant(const T& m) {
  check_square("determinant", "m", m);
  if (m.size() == 0) {
    return 1;
  }

  double det_val = m.val().determinant();
  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arena_m = m;
  arena_matrix<Eigen::MatrixXd> arena_m_inv_t = m.val().inverse().transpose();

  var det = det_val;

  reverse_pass_callback([=]() mutable {
    arena_m.adj() += (det.adj() * det_val) * arena_m_inv_t;
  });

  return det;
}

}  // namespace math
}  // namespace stan
#endif
