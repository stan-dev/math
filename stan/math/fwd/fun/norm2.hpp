#ifndef STAN_MATH_FWD_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_FWD_FUN_LOG_SUM_EXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/norm2.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Compute the L2 norm of the specified vector of values.
 *
 * @tparam T Type of input vector.
 * @param[in] x Vector of specified values.
 * @return L2 norm of x.
 */
template <typename Container,
          require_container_st<is_fvar, Container>* = nullptr>
inline auto norm2(const Container& x) {
  return apply_vector_unary<ref_type_t<T>>::reduce(
      to_ref(x), [&](const auto& v) {
        using T_fvar_inner = typename value_type_t<decltype(v)>::Scalar;
        T_fvar_inner res = norm2(v.val().array());
        return fvar<T_fvar_inner>(res, v.d().array() * (v.val().array() / res));
      });
}

}  // namespace math
}  // namespace stan
#endif
