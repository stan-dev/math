#ifndef STAN_MATH_REV_META_REV_MATRIX_TYPE_HPP
#define STAN_MATH_REV_META_REV_MATRIX_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta/is_var.hpp>

namespace stan {

/**
 * Determines a return type for a function that accepts given inputs and wants
 * to return a matrix (or vector or row vector) with given compile time number
 * of rows and columns.
 *
 * Returns `var_value` if any of the inputs are `var_value`, otherwise an Eigen
 * matrix of `var` or `double`.
 * @tparam Rows number of compile time rows
 * @tparam Cols number of compile time columns
 * @tparam Inputs types of arguments the function accepts
 */
template <int Rows, int Cols, typename... Inputs>
struct rev_matrix_type {
  using type = std::conditional_t<
      math::disjunction<math::conjunction<
          is_var<Inputs>, is_eigen<value_type_t<Inputs>>>...>::value,
      math::var_value<Eigen::Matrix<double, Rows, Cols>>,
      Eigen::Matrix<return_type_t<Inputs...>, Rows, Cols>>;
};

template <int Rows, int Cols, typename... Inputs>
using rev_matrix_t = typename rev_matrix_type<Rows, Cols, Inputs...>::type;

}  // namespace stan

#endif
