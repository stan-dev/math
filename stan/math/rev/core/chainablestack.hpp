#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
#include <tuple>
namespace stan {
namespace math {

template <typename T, typename>
class vari_value;
class vari_base;
class chainable_alloc;
using ChainableStack
    = AutodiffStackSingleton<vari_base, chainable_alloc>;

}  // namespace math
}  // namespace stan
#endif
