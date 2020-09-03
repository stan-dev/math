#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
namespace stan {
namespace math {
class chainable_alloc;
class vari_base;
using ChainableStack = AutodiffStackSingleton<vari_base, chainable_alloc>;

}  // namespace math
}  // namespace stan
#endif
