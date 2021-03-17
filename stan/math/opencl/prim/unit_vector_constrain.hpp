#ifndef STAN_MATH_OPENCL_UNIT_VECTOR_CONSTRAIN_BLOCK_HPP
#define STAN_MATH_OPENCL_UNIT_VECTOR_CONSTRAIN_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the given free vector.
 *
 * See <a
 * href="https://en.wikipedia.org/wiki/N-sphere#Generating_random_points">the
 * Wikipedia page on generating random points on an N-sphere</a>.
 *
 * @tparam T_x type of the vector
 * @param x vector of unrestricted variables
 * @return Unit length vector
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline auto unit_vector_constrain(T_x&& x) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "x", x);
  value_type_t<T_x> SN = dot_self(x);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  return elt_divide(x, sqrt(SN));
}

/**
 * Return the unit length vector corresponding to the given free vector.
 *
 * See <a
 * href="https://en.wikipedia.org/wiki/N-sphere#Generating_random_points">the
 * Wikipedia page on generating random points on an N-sphere</a>.
 *
 * @tparam T_x type of the vector
 * @param x vector of unrestricted variables
 * @param[in,out] lp log probability to increment
 * @return Unit length vector
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline auto unit_vector_constrain(T_x&& x, double& lp) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "x", x);
  value_type_t<T_x> SN = dot_self(x);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  lp -= 0.5 * SN;
  return elt_divide(x, sqrt(SN));
}

}  // namespace math
}  // namespace stan
#endif
#endif
