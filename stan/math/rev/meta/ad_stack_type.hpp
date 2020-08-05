#ifndef AD_STACK_TYPE_HPP
#define AD_STACK_TYPE_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/rev/core/ad_allocator.hpp>

namespace stan {
namespace math {
//forward declaration
template<typename MatrixType>
class AD_stack_matrix;
}

namespace internal{
template<typename T, typename=void>
struct AD_stack_type_impl{
    using type = T;
};

template<typename T>
struct AD_stack_type_impl<std::vector<T>>{
    using T_ad = typename AD_stack_type_impl<std::decay_t<T>>::type;
    using type = std::vector<T_ad, math::AD_allocator<T_ad>>;
};

template<typename T>
struct AD_stack_type_impl<T, require_eigen_t<T>>{
    using type = math::AD_stack_matrix<plain_type_t<T>>;
};
}

/**
 * Determines a type that can be used in place of `T` that does any dynamic allocations on the AD stack. This way resulting types are trivially destructible and can be used in vari classes.
 * (only works for POD types, `std::vector`s and Eigen types)
 */
template <typename T>
using AD_stack_t = typename internal::AD_stack_type_impl<std::decay_t<T>>::type;

}

#endif // AD_STACK_TYPE_HPP
