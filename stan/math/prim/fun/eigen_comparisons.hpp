#ifndef STAN_MATH_PRIM_FUN_EIGEN_COMPARISONS_HPP
#define STAN_MATH_PRIM_FUN_EIGEN_COMPARISONS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Add support for comparisons involving Eigen types with different scalars,
 * where one of the scalars is an autodiff type. This includes comparisons of
 * an Eigen type and a scalar.
 * @param OPERATOR name of the operator function to implement
 * @param OP operator to use for comparison of values
 **/
#define ADD_MIXED_AUTODIFF_SCALAR_COMPARISON(OPERATOR, OP) \
  template <typename T_a, typename T_b,                    \
            require_any_eigen_t<T_a, T_b>* = nullptr,      \
            require_any_st_autodiff<T_a, T_b>* = nullptr,  \
            require_not_st_same<T_a, T_b>* = nullptr>      \
  auto OPERATOR(const T_a& a, const T_b& b) {              \
    return value_of(a) OP value_of(b);                     \
  }

ADD_MIXED_AUTODIFF_SCALAR_COMPARISON(operator<, <);
ADD_MIXED_AUTODIFF_SCALAR_COMPARISON(operator<=, <=);
ADD_MIXED_AUTODIFF_SCALAR_COMPARISON(operator>, >);
ADD_MIXED_AUTODIFF_SCALAR_COMPARISON(operator>=, >=);
ADD_MIXED_AUTODIFF_SCALAR_COMPARISON(operator==, ==);
ADD_MIXED_AUTODIFF_SCALAR_COMPARISON(operator!=, !=);

#undef ADD_MIXED_AUTODIFF_SCALAR_COMPARISON

}  // namespace math
}  // namespace stan

#endif
