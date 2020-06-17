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
using zeroing_stacks = std::tuple<std::vector<vari_value<double>*>,
  std::vector<vari_value<float>*>,
  std::vector<vari_value<long double>*>>;
using ChainableStack = AutodiffStackSingleton<vari_base, chainable_alloc, zeroing_stacks>;

}  // namespace math
}  // namespace stan
#endif
