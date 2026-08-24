#ifndef STAN_MATH_PRIM_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_PRIM_FUN_LOG_SOFTMAX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the softmax of the specified vector or 
 * matrix, or of each vector or matrix in a container. For a matrix, the 
 * log-softmax is taken over all elements.
 *  *
 * \f$
 * \log \mbox{softmax}(y)
 * \ = \ y - \log \sum_{k=1}^K \exp(y_k)
 * \ = \ y - \mbox{log\_sum\_exp}(y).
 * \f$
 *
 * For the log softmax function, the entries in the Jacobian are
 * \f$
 * \frac{\partial}{\partial y_m} \log\mbox{softmax}(y)[k]
 * = \left\{
 * \begin{array}{ll}
 * 1 - \mbox{softmax}(y)[m]
 * & \mbox{ if } m = k, \mbox{ and}
 * \\[6pt]
 * -\mbox{softmax}(y)[m]
 * & \mbox{ if } m \neq k.
 * \end{array}
 * \right.
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
inline auto log_softmax(Container&& x) {
  return make_holder(
      [](auto&& a) {
        return apply_vector_unary<ref_type_t<Container>>::apply(
            std::forward<decltype(a)>(a),
            [](auto&& v) -> plain_type_t<decltype(v)> {
              if (v.size() == 0) {
                return v;
              }
              return (v.array() - log_sum_exp(v)).matrix();
            });
      },
      to_ref(std::forward<Container>(x)));
}

}  // namespace math
}  // namespace stan
#endif
