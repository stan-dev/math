#ifndef STAN_MATH_FWD_FUN_NORM2_HPP
#define STAN_MATH_FWD_FUN_NORM2_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/norm2.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Compute the L2 norm of the specified vector of values.
 *
 * @tparam T Type of input vector.
 * @param[in] x Vector of specified values.
 * @return L2 norm of x.
 */
template <typename Container, require_eigen_vt<is_fvar, Container>* = nullptr>
inline auto norm2(const Container& x) {
  return apply_vector_unary<ref_type_t<Container>>::reduce(
      to_ref(x), [&](const auto& v) {
        using T_fvar_inner = typename value_type_t<decltype(v)>::Scalar;
        T_fvar_inner res = norm2(v.val());
        return fvar<T_fvar_inner>(res,
                                  v.d().cwiseProduct((v.val() / res)).sum());
      });
}

}  // namespace math
}  // namespace stan
#endif
