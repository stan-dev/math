#ifndef STAN_MATH_FWD_FUN_NORM1_HPP
#define STAN_MATH_FWD_FUN_NORM1_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/norm1.hpp>

namespace stan {
namespace math {

/**
 * Compute the L1 norm of the specified vector of values.
 *
 * @tparam T Type of input vector.
 * @param[in] x Vector of specified values.
 * @return L1 norm of x.
 */
template <typename Container, require_eigen_vt<is_fvar, Container>* = nullptr>
inline auto norm1(const Container& x) {
  return apply_vector_unary<ref_type_t<Container>>::reduce(
      to_ref(x), [&](const auto& v) {
        using T_fvar_inner = typename value_type_t<decltype(v)>::Scalar;
        return fvar<T_fvar_inner>(norm1(v.val()),
                                  v.d().cwiseProduct(sign(v.val())).sum());
      });
}

}  // namespace math
}  // namespace stan
#endif
