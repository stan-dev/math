#ifndef STAN_MATH_PRIM_FUN_SOFTMAX_HPP
#define STAN_MATH_PRIM_FUN_SOFTMAX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified vector or matrix, or of each
 * vector or matrix in a container. For a matrix, the softmax is
 * taken over all elements.
 *
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
 * @tparam Container type of input: an Eigen vector, Eigen matrix, 
 *   `std::vector` of vectors or matrices, or nested container whose scalar 
 *   type is arithmetic
 * @param x vector, matrix, or container to transform.
 * @return softmax of the input, preserving the container structure; an empty
 *   result if any input vector or matrix is empty.
 */
template <typename Container, require_st_arithmetic<Container>* = nullptr,
          require_container_t<Container>* = nullptr>
inline auto softmax(Container&& x) {
  return make_holder(
      [](auto&& a) {
        return apply_vector_unary<ref_type_t<Container>>::apply(
            std::forward<decltype(a)>(a),
            [](auto&& v) -> plain_type_t<decltype(v)> {
              if (v.size() == 0) {
                return v;
              }
              const auto theta = (v.array() - v.maxCoeff()).exp();
              return (theta / theta.sum()).matrix();
            });
      },
      to_ref(std::forward<Container>(x)));
}

}  // namespace math
}  // namespace stan
#endif
