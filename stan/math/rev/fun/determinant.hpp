#ifndef STAN_MATH_REV_FUN_DETERMINANT_HPP
#define STAN_MATH_REV_FUN_DETERMINANT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline var determinant(const T& m) {
  check_square("determinant", "m", m);

  if (m.size() == 0) {
    return 1;
  }

  const auto& m_ref = to_ref(m);
  const auto& m_val = value_of(m_ref);

  double det_val = m_val.determinant();
  arena_t<plain_type_t<T>> arena_m = m_ref;
  arena_t<Eigen::MatrixXd> arena_m_inv_t = m_val.inverse().transpose();

  var det = det_val;

  reverse_pass_callback([arena_m, det, arena_m_inv_t]() mutable {
    arena_m.adj() += (det.adj() * det.val()) * arena_m_inv_t;
  });

  return det;
}

}  // namespace math
}  // namespace stan
#endif
