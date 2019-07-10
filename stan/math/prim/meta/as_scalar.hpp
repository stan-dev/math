#ifndef STAN_MATH_PRIM_META_AS_SCALAR_HPP
#define STAN_MATH_PRIM_META_AS_SCALAR_HPP

#include <Eigen/Dense>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Converts input to a scalar. As this is not possible for matrices, arrays or
 * Eigen expressions it always throws. This is intended to never be called, only
 * used in templated functions in branches that will be optimized out - to
 * prevent compiler from complaining about expressions with incompatible types.
 * @tparam Derived Type of input Eigen expression.
 * @param a Input expression
 * @throws runtime_error Always throws
 * @return Never returns
 */
template <typename Derived>
inline double as_scalar(const Eigen::DenseBase<Derived>& a) {
  throw std::runtime_error("A matrix can not be used as a scalar!");
}

/**
 * Converts input to a scalar. As this is not possible for vectors it always
 * throws. This is intended to never be called, only used in templated functions
 * in branches that will be optimized out - to prevent compiler from complaining
 * about expressions with incompatible types.
 * @param a Input expression
 * @throws runtime_error Always throws
 * @return Never returns
 */
template <typename T>
inline double as_scalar(const std::vector<T>& a) {
  throw std::runtime_error("A vector can not be used as a scalar!");
}
  
/**
 * Converts input to a scalar. For scalar arguments this is an identity
 * function.
 * @param a Input value
 * @return Same value
 */
inline double as_scalar(double a) { return a; }

/**
 * Converts input to a scalar. For scalar arguments this is an identity
 * function.
 * @param a Input value
 * @return Same value
 */
inline int as_scalar(int a) { return a; }
  
}  // namespace math
}  // namespace stan

#endif
