#ifndef STAN_MATH_REV_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_REV_FUN_LOG_SOFTMAX_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the log softmax of the specified vector or container of vectors.
 *
 * The gradient calculations are unfolded.
 *
 * @tparam T Type of input vector or matrix.
 * @param[in] x Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */
template <typename T, require_container_st<is_var, T>* = nullptr>
inline auto log_softmax(const T& x) {
  return apply_vector_unary<ref_type_t<T>>::apply(
      to_ref(x), [](const auto& alpha) {
        using ret_type = plain_type_t<decltype(alpha)>;

        check_nonzero_size("log_softmax", "alpha", alpha);

        auto arena_alpha = to_arena(to_ref(alpha));
        const auto& alpha_val = to_ref(value_of(arena_alpha).array());
        const auto& theta = to_ref(alpha_val - alpha_val.maxCoeff());

        arena_t<ret_type> res = theta.array() - log(theta.exp().sum());

        reverse_pass_callback([arena_alpha, res]() mutable {
          arena_alpha.adj()
              += res.adj()
                 - (res.adj().sum() * res.val().array().exp()).matrix();
        });

        return ret_type(res);
      });
}

}  // namespace math
}  // namespace stan
#endif
