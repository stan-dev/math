#ifndef STAN_MATH_PRIM_SCAL_FUN_INV_SQUARE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_INV_SQUARE_HPP

#include <stan/math/prim/scal/fun/inv.hpp>
#include <stan/math/prim/scal/fun/square.hpp>

namespace stan {
namespace math {

inline double inv_square(double x) { return inv(square(x)); }
}  // namespace math
}  // namespace stan

#endif
