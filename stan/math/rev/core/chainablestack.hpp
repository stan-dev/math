#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
#include <stan/math/rev/core/make_vari_variadic.hpp>
#include <boost/variant2/variant.hpp>
namespace stan {
namespace math {
class chainable_alloc;
using ChainableStack = AutodiffStackSingleton<vari_variant, chainable_alloc>;

}  // namespace math
}  // namespace stan
#endif
