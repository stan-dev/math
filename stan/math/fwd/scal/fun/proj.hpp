#ifndef STAN_MATH_FWD_SCAL_FUN_PROJ_HPP
#define STAN_MATH_FWD_SCAL_FUN_PROJ_HPP

#include <stan/math/fwd/cplx.hpp>
#include <stan/math/fwd/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/proj.hpp>

namespace std {

/**
 * Return a projection of a complex number onto
 * the Riemann sphere.
 *
 * This overloads an erroneous definition in
 * libstdc++ for general types.
 *
 * @tparam T auto diff variable type
 * @param[in] t Variable input.
 * @return projection onto Riemann sphere
 */
template <class T>
inline std::complex<stan::math::fvar<T>> proj(
    std::complex<stan::math::fvar<T>> const& t) {
  return stan::math::proj(t);
}

}  // namespace std
#endif
