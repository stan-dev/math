#ifndef STAN_MATH_PRIM_FUN_COS_HPP
#define STAN_MATH_PRIM_FUN_COS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/cosh.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/i_times.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap `cos()` so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x angle in radians
 * @return Cosine of x.
 */
struct cos_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    using std::cos;
    return cos(x);
  }
};

/**
 * Returns the elementwise `cos()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x angles in radians
 * @return Cosine of each value in x.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto cos(const Container& x) {
  return apply_scalar_unary<cos_fun, Container>::apply(x);
}

/**
 * Version of `cos()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Cosine of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto cos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [&](const auto& v) { return v.array().cos(); });
}

namespace internal {
/**
 * Return the cosine of the complex argument.
 *
 * @tparam T value type of argument
 * @param[in] z argument
 * @return cosine of the argument
 */
template <typename T>
inline std::complex<T> complex_cos(const std::complex<T>& z) {
  return cosh(i_times(z));
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
