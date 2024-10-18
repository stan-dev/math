
#ifndef STAN_MATH_VECTORIZED_CONSTRAINT_CHOLESKY_FACTOR_CONSTRAIN_HPP
#define STAN_MATH_VECTORIZED_CONSTRAINT_CHOLESKY_FACTOR_CONSTRAIN_HPP
#include <stan/math/prim/constraint/cholesky_factor_constrain.hpp>
namespace stan {
namespace math {

/**
 * Return the Cholesky factor of the specified size read from the specified
 * vector. A total of (N choose 2) + N + N * (M - N) free parameters are
 * required to read an M by N Cholesky factor. If the `Jacobian` parameter is
 * `true`, the log density accumulator is incremented with the log absolute
 * Jacobian determinant of the transform.  All of the transforms are specified
 * with their Jacobians in the *Stan Reference Manual* chapter Constraint
 * Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param x Vector of unconstrained values
 * @param M number of rows
 * @param N number of columns
 * @param[in,out] lp log density accumulator
 * @return Cholesky factor
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto cholesky_factor_constrain(const T& x, int M, int N,
                                      return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(x, [&lp, M, N](auto&& v) {
    return cholesky_factor_constrain<Jacobian>(v, M, N, lp);
  });
}


} // namespace math
} // namespace stan
#endif 

