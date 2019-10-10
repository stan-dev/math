#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

using ChainableStack = AutodiffStackSingleton<vari, chainable_alloc>;

}  // namespace math
}  // namespace stan
#endif
