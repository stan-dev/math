#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

template <typename T>
class vari_type;
class chainable_alloc;

using ChainableStack = AutodiffStackSingleton<vari_type<double>, chainable_alloc>;

}  // namespace math
}  // namespace stan
#endif
