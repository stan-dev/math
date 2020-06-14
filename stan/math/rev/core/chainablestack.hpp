#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
#include <boost/variant2/variant.hpp>
namespace stan {
namespace math {

template <typename T, typename = void>
class vari_value;
class vari_base;
class chainable_alloc;

using vari_variant
    = boost::variant2::variant<vari_value<double>*, vari_value<float>*, vari_value<long double>*>;

using ChainableStack = AutodiffStackSingleton<vari_variant, chainable_alloc>;

}  // namespace math
}  // namespace stan
#endif
