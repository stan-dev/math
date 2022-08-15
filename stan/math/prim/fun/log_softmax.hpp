#ifndef STAN_MATH_PRIM_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_PRIM_FUN_LOG_SOFTMAX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/softmax.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the softmax of the specified
 * vector.
 *
 * \f$
 * \log \mbox{softmax}(y)
 * \ = \ y - \log \sum_{k=1}^K \exp(y_k)
 * \ = \ y - \mbox{log\_sum\_exp}(y).
 * \f$
 *
 * @tparam Container type of input vector to transform
 * @param[in] x vector to transform
 * @return log unit simplex result of the softmax transform of the vector.
 */
template <typename V,
          require_eigen_col_vector_vt<std::is_arithmetic, V>* = nullptr>
inline auto log_softmax(const V& v) {
  using stan::math::log;
  using stan::math::softmax;
  check_nonzero_size("log_softmax", "v", v);
  auto z = softmax(v);
  auto y = log(z);
  return y.eval();
}

}  // namespace math
}  // namespace stan
#endif
