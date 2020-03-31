#ifndef STAN_MATH_PRIM_FUN_SOFTMAX_HPP
#define STAN_MATH_PRIM_FUN_SOFTMAX_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified Eigen expression/object,
 * std::vector, or container of these.
 *
 * <p>
 * \f$
 * \mbox{softmax}(y)
 * = \frac{\exp(y)}
 * {\sum_{k=1}^K \exp(y_k)},
 * \f$
 *
 * <p>The entries in the Jacobian of the softmax function are given by
 * \f$
 * \begin{array}{l}
 * \displaystyle
 * \frac{\partial}{\partial y_m} \mbox{softmax}(y)[k]
 * \\[8pt]
 * \displaystyle
 * \mbox{ } \ \ \ = \left\{
 * \begin{array}{ll}
 * \mbox{softmax}(y)[k] \times (1 - \mbox{softmax}(y)[m])
 * & \mbox{ if } m = k, \mbox{ and}
 * \\[6pt]
 * -\mbox{softmax}(y)[k] \times \mbox{softmax}(y)[m]
 * & \mbox{ if } m \neq k.
 * \end{array}
 * \right.
 * \end{array}
 * \f$
 *
 * @tparam Container type of container
 * @param[in] x Container (or container of containers) to transform.
 * @return Unit simplex result of the softmax transform of the container.
 */
template <typename Container,
          require_vector_st<std::is_arithmetic, Container>...>
inline auto softmax(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [](const auto& v) {
    const Eigen::Ref<const plain_type_t<decltype(v)>>& v_ref = v;
    check_nonzero_size("softmax", "v", v_ref);
    plain_type_t<decltype(v)> theta(v_ref.size());
    theta = (v_ref.array() - v_ref.maxCoeff()).exp();
    return (theta.array() / theta.sum()).eval();
  });
}

}  // namespace math
}  // namespace stan

#endif
