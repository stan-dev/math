#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_COMMON_RETURN_SCALAR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_COMMON_RETURN_SCALAR_HPP
#include <type_traits>

namespace stan {
namespace math {
/**
 * @ingroup type_traits
 * Wrapper for `std::common_type_t`
 */
template <typename... Types>
using common_scalar_t = typename std::common_type_t<
    typename std::remove_reference_t<Types>::Scalar...>;
}  // namespace math
}  // namespace stan

#endif
