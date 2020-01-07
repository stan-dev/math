#ifndef STAN_MATH_REV_FUN_STAN_PRINT_HPP
#define STAN_MATH_REV_FUN_STAN_PRINT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <ostream>

namespace stan {
namespace math {

inline void stan_print(std::ostream* o, const var& x) { *o << x.val(); }

}  // namespace math
}  // namespace stan
#endif
