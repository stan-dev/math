#ifndef STAN_MATH_REV_CORE_SET_ZERO_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_SET_ZERO_ADJOINTS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/functor/iter_tuple_nested.hpp>
#include <type_traits>
namespace stan::math::internal {
/**
 * Set all adjoints of the output to zero.
 * @param[in, out] output The output whose adjoints will be set to zero
 */
template <typename Output>
inline void set_zero_adjoint(Output&& output) {
  if constexpr (is_all_arithmetic_scalar_v<Output>) {
    return;
  } else {
    return iter_tuple_nested(
        [](auto&& output_i) {
          using output_i_t = std::decay_t<decltype(output_i)>;
          if constexpr (is_all_arithmetic_scalar_v<output_i_t>) {
            return;
          } else if constexpr (is_std_vector<output_i_t>::value) {
            for (Eigen::Index i = 0; i < output_i.size(); ++i) {
              output_i[i].adj() = 0;
            }
          } else if constexpr (is_eigen_v<output_i_t>) {
            output_i.adj().setZero();
          } else if constexpr (is_stan_scalar_v<output_i_t>) {
            output_i.adj() = 0;
          } else {
            static_assert(
                sizeof(std::decay_t<output_i_t>*) == 0,
                "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
                "not able to deduce the actions needed for the given type. "
                "This is an internal error, please report it: "
                "https://github.com/stan-dev/math/issues");
          }
        },
        std::forward<Output>(output));
  }
}
}  // namespace stan::math::internal

#endif
