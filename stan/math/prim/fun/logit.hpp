#ifndef STAN_MATH_PRIM_FUN_LOGIT_HPP
#define STAN_MATH_PRIM_FUN_LOGIT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the log odds of the argument.
 *
 * The logit function is defined as for \f$x \in [0, 1]\f$ by
 * returning the log odds of \f$x\f$ treated as a probability,
 *
 * \f$\mbox{logit}(x) = \log \left( \frac{x}{1 - x} \right)\f$.
 *
 * The inverse to this function is <code>inv_logit</code>.
 *
 *
 \f[
 \mbox{logit}(x) =
 \begin{cases}
 \textrm{NaN}& \mbox{if } x < 0 \textrm{ or } x > 1\\
 \ln\frac{x}{1-x} & \mbox{if } 0\leq x \leq 1 \\[6pt]
 \textrm{NaN} & \mbox{if } x = \textrm{NaN}
 \end{cases}
 \f]

 \f[
 \frac{\partial\, \mbox{logit}(x)}{\partial x} =
 \begin{cases}
 \textrm{NaN}& \mbox{if } x < 0 \textrm{ or } x > 1\\
 \frac{1}{x-x^2}& \mbox{if } 0\leq x\leq 1 \\[6pt]
 \textrm{NaN} & \mbox{if } x = \textrm{NaN}
 \end{cases}
 \f]
 *
 * @param u argument
 * @return log odds of argument
 */
inline double logit(double u) {
  using std::log;
  return log(u / (1 - u));
}

/**
 * Return the log odds of the argument.
 *
 * @param u argument
 * @return log odds of argument
 */
inline double logit(int u) { return logit(static_cast<double>(u)); }

/**
 * Structure to wrap logit() so it can be vectorized.
 */
struct logit_fun {
  /**
   * Return the log odds of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return log odds of the argument
   */
  template <typename T>
  static inline T fun(const T& x) {
    return logit(x);
  }
};

/**
 * Return the elementwise application of <code>logit()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam Container type of container
 * @param x container
 * @return elementwise logit of container elements
 */
template <
    typename Container,
    require_not_container_st<std::is_arithmetic, Container>* = nullptr,
    require_not_var_matrix_t<Container>* = nullptr,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<Container>* = nullptr>
inline auto logit(const Container& x) {
  return apply_scalar_unary<logit_fun, Container>::apply(x);
}

/**
 * Version of logit() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return the logit of each variable in the container.
 *
 * Note: The return must be evaluated otherwise the Ref object falls out
 * of scope
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto logit(const Container& x) {
  return make_holder(
      [](const auto& v_ref) {
        return apply_vector_unary<ref_type_t<Container>>::apply(
            v_ref,
            [](const auto& v) { return (v.array() / (1 - v.array())).log(); });
      },
      to_ref(x));
}

}  // namespace math
}  // namespace stan

#endif
