#ifndef STAN_MATH_REV_SCAL_FUN_PROJ_HPP
#define STAN_MATH_REV_SCAL_FUN_PROJ_HPP

#include <stan/math/rev/cplx.hpp>
#include <stan/math/rev/scal/fun/is_inf.hpp> 
#include <stan/math/prim/scal/fun/proj.hpp>

namespace std {

/**
 * Return a projection of a complex number onto
 * the Riemann sphere.
 *
 * This overloads an erroneous definition in
 * libstdc++ for general types.
 *
 * @param[in] t Variable input.
 * @return projection onto Riemann sphere
 */

inline std::complex<stan::math::var>
proj(std::complex<stan::math::var> const& t) {
  return stan::math::proj(t);
}

}  // namespace std
#endif
