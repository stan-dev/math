#ifndef STAN_MATH_REV_META_IS_VARI_HPP
#define STAN_MATH_REV_META_IS_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier and pointer is a vari.
 */
template <typename T>
struct is_vari<T, std::enable_if_t<
  std::is_base_of<math::vari_base,
    std::remove_pointer_t<std::decay_t<T>>>::value>> : std::true_type {};

}  // namespace stan
#endif
