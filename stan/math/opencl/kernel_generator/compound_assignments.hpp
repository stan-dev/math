#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_COMPOUND_ASSIGNMENTS_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_COMPOUND_ASSIGNMENTS_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Adds a compound assignment operator for kernel generator expressions.
 * @param independent independent (non-compound) operator
 * @param compound compound operator
 */
#define ADD_KG_COMPOUND_ASSIGMENT_OPERATOR(independent, compound) \
  template <typename T1, typename T2,                             \
            require_kernel_expression_lhs_t<T1>* = nullptr,       \
            require_all_kernel_expressions_t<T2>* = nullptr>      \
  T1 operator compound(T1&& a, T2&& b) {                          \
    return a = a independent b;                                   \
  }

ADD_KG_COMPOUND_ASSIGMENT_OPERATOR(+, +=)
ADD_KG_COMPOUND_ASSIGMENT_OPERATOR(-, -=)
ADD_KG_COMPOUND_ASSIGMENT_OPERATOR(*, *=)

#undef ADD_KG_COMPOUND_ASSIGMENT_OPERATOR

}  // namespace math
}  // namespace stan

#endif  // COMPOUND_ASSIGNMENTS_HPP
