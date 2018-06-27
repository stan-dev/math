#ifndef STAN_MATH_PRIM_SCAL_FUN_POLAR_HPP
#define STAN_MATH_PRIM_SCAL_FUN_POLAR_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/rm_zeroing.hpp>
#include <complex>
#include <cmath>
#include <type_traits>

namespace std {
/**
 * Returns a complex number with magnitude r
 * and phase theta
 *
 * This is in the std namespace because the 
 * relevant call we are catching is fully
 * qualified with std, making ADL useless.
 *
 * @tparam R magnitude type
 * @tparam TH phase type
 * @param r magnitude
 * @param theta angle, radians
 * @return complex with magnitude r & phase theta
 */
template <class R, class TH, class RT = typename
 stan::return_type<stan::rm_zeroing_t<R>,
  stan::rm_zeroing_t<TH>>::type,
 std::enable_if_t<
 !std::is_same<R, TH>::value>* =nullptr>
inline std::complex<RT>
polar(R const& r, TH const& theta) {
 return polar(RT(r),RT(theta));
}

}  // namespace std
#endif
