#ifndef STAN_MATH_PRIM_FUN_ATAN_HPP
#define STAN_MATH_PRIM_FUN_ATAN_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/atanh.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/i_times.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap \c atan() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Arctan of x in radians.
 */
struct atan_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    using std::atan;
    return atan(x);
  }
};

/**
 * Returns the elementwise \c atan() of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return Arctan of each value in x, in radians.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto atan(const Container& x) {
  return apply_scalar_unary<atan_fun, Container>::apply(x);
}

/**
 * Version of atan() that accepts std::vectors, Eigen Matrix/Array objects,
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Elementwise atan of members of container.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto atan(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().atan(); });
}

namespace internal {
/**
 * Return the arc tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return arc tangent of the argument
 */
template <typename V>
inline std::complex<V> complex_atan(const std::complex<V>& z) {
  return neg_i_times(atanh(i_times(z)));
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
