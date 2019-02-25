#ifndef STAN_MATH_PRIM_FUN_MODULUS_HPP
#define STAN_MATH_PRIM_FUN_MODULUS_HPP

#include <stanh/prim/err/domain_error.hpp>
#include <stanh/prim/meta/likely.hpp>
#include <cstddef>
#include <cstdlib>

namespace stan {
namespace math {

inline int modulus(int x, int y) {
  if (unlikely(y == 0))
    domain_error("modulus", "divisor is", 0, "");
  return x % y;
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_MODULUS_HPP
#define STAN_MATH_PRIM_FUN_MODULUS_HPP

#include <stanh/prim/err/domain_error.hpp>
#include <stanh/prim/meta/likely.hpp>
#include <cstddef>
#include <cstdlib>

namespace stan {
namespace math {

inline int modulus(int x, int y) {
  if (unlikely(y == 0))
    domain_error("modulus", "divisor is", 0, "");
  return x % y;
}

}  // namespace math
}  // namespace stan
#endif
