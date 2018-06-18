#ifndef STAN_MATH_PRIM_SCAL_FUN_COMPLEX_PROMOTE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_COMPLEX_PROMOTE_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/rm_complex.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Complex auto-diff (AD) promotion for arguments
 * to binary operators. Handles cases where types
 * may or may not be wrapped in complex.
 *
 * @tparam N non-deduced input type
 * @tparam D deduced input type
 * @param d complex argument
 * @return value of promoted type
 */
template <class N, class D, class P = typename 
 return_type<rm_complex_t<N>,
             rm_complex_t<D>>::type>
inline auto
complex_promote(std::complex<D> const& d) {
  return std::complex<P>(d.real(), d.imag());
}

/**
 * Complex auto-diff (AD) promotion for arguments
 * to binary operators. Handles cases where types
 * may or may not be wrapped in complex.
 *
 * @tparam N non-deduced input type
 * @tparam D deduced input type
 * @param d argument
 * @return value of promoted type
 */
template <class N, class D, class P = typename 
 return_type<rm_complex_t<N>,
             rm_complex_t<D>>::type>
inline auto
complex_promote(D const& d) {
  return std::complex<P>(d);
}

}  // namespace math
}  // namespace stan
#endif
