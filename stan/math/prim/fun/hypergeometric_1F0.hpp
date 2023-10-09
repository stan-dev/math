#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_1F0_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_1F0_HPP

#include <stan/math/prim/meta.hpp>
#include <boost/math/special_functions/hypergeometric_1F0.hpp>

namespace stan {
namespace math {

template <typename Ta, typename Tz,
          require_all_arithmetic_t<Ta, Tz>* = nullptr>
auto hypergeometric_1f0(const Ta& a, const Tz& z) {
  return boost::math::hypergeometric_1F0(a, z);
}

}  // namespace math
}  // namespace stan
#endif
