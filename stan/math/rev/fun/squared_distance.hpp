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
  check_finite("squared_distance", "a", a);
  check_finite("squared_distance", "b", b);

  double diff = value_of(a) - value_of(b);
  var res = diff * diff;

  reverse_pass_callback([a, b, res]() mutable {
    double diff = value_of(a) - value_of(b);
    if(!is_constant<T1>::value)
      forward_as<var>(a).adj() += 2.0 * res.adj() * diff;
    if(!is_constant<T2>::value)
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

  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);

  double res_val = 0.0;
  for(size_t i = 0; i < A.size(); ++i) {
    double diff = value_of(A_ref[i]) - value_of(B_ref[i]);
    res_val += diff * diff;
  }

  var res = res_val;

  auto arena_A = to_arena(A_ref);
  auto arena_B = to_arena(B_ref);

  reverse_pass_callback([arena_A, arena_B, res]() mutable {
    for(size_t i = 0; i < arena_A.size(); ++i) {
      double diff = value_of(arena_A[i]) -
	value_of(arena_B[i]);
      if(!is_constant<T1>::value)
	forward_as<var>(arena_A[i]).adj() += 2.0 * res.adj() * diff;
      if(!is_constant<T2>::value)
	forward_as<var>(arena_B[i]).adj() -= 2.0 * res.adj() * diff;
    }
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
