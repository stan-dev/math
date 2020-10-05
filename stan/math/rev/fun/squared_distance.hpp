#ifndef STAN_MATH_REV_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_REV_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the squared distance.
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr,
          require_any_var_t<T1, T2>* = nullptr>
inline var squared_distance(const T1& a, const T2& b) {
  double diff = value_of(a) - value_of(b);
  var res = squared_distance(value_of(a), value_of(b));
  reverse_pass_callback([a, b, diff, res]() mutable {
    if (!is_constant<T1>::value)
      forward_as<var>(a).adj() += 2.0 * res.adj() * diff;

    if (!is_constant<T2>::value)
      forward_as<var>(b).adj() += -2.0 * res.adj() * diff;
  });
  return res;
}

template <typename T1, typename T2,
          require_all_eigen_vector_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr>
inline var squared_distance(const T1& A, const T2& B) {
  check_matching_sizes("squared_distance", "A", A, "B", B);

  if (A.size() == 0)
    return 0.0;

  using A_ref_t = ref_type_t<T1>;
  using B_ref_t = ref_type_t<T2>;

  A_ref_t A_ref = A;
  B_ref_t B_ref = B;

  auto A_col = as_column_vector_or_scalar(A_ref);
  auto B_col = as_column_vector_or_scalar(B_ref);

  auto arena_diff_val = to_arena(value_of(A_col) - value_of(B_col));

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> arena_A;
  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> arena_B;

  if (!is_constant<T1>::value) {
    arena_A = A_col;
  }

  if (!is_constant<T2>::value) {
    arena_B = B_col;
  }

  var res = arena_diff_val.squaredNorm();

  reverse_pass_callback([arena_A, arena_B, arena_diff_val, res]() mutable {
    if (!is_constant<T1>::value)
      arena_A.adj() += 2.0 * res.adj() * arena_diff_val;

    if (!is_constant<T2>::value)
      arena_B.adj() += -2.0 * res.adj() * arena_diff_val;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
