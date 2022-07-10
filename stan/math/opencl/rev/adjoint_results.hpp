#ifndef STAN_MATH_OPENCL_REV_ADJOINT_RESULTS_HPP
#define STAN_MATH_OPENCL_REV_ADJOINT_RESULTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/multi_result_kernel.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents results that are adjoints of vars in kernel generrator
 * expressions.
 */
template <typename... T_results>
class adjoint_results_cl : protected results_cl<T_results...> {
 public:
  /**
   * Constructor.
   * @param results results that will be calculated in same kernel
   */
  explicit adjoint_results_cl(T_results&&... results)
      : results_cl<T_results...>(std::forward<T_results>(results)...) {}

  /**
   * Incrementing \c adjoint_results_cl object by \c expressions_cl object
   * executes one or two kernels that evaluate expressions and increments
   * results' adjoints by those expressions. Non-var results are ignored.
   * Scalar results's adjoints get incremented by the sum of respective
   * expressions.
   * @tparam T_expressions types of expressions
   * @param exprs expressions
   */
  template <typename... T_expressions,
            typename = std::enable_if_t<sizeof...(T_results)
                                        == sizeof...(T_expressions)>>
  void operator+=(const expressions_cl<T_expressions...>& exprs) {
    index_apply<sizeof...(T_expressions)>([&](auto... Is) {
      auto scalars = std::tuple_cat(select_scalar_assignments(
          std::get<Is>(this->results_), std::get<Is>(exprs.expressions_))...);
      auto nonscalars_tmp = std::tuple_cat(select_nonscalar_assignments(
          std::get<Is>(this->results_), std::get<Is>(exprs.expressions_))...);

      index_apply<std::tuple_size<decltype(nonscalars_tmp)>::value>(
          [&](auto... Is_nonscal) {
            auto nonscalars = std::make_tuple(std::make_pair(
                std::get<Is_nonscal>(nonscalars_tmp).first,
                std::get<Is_nonscal>(nonscalars_tmp).first
                    + std::get<Is_nonscal>(nonscalars_tmp).second)...);

            index_apply<std::tuple_size<decltype(scalars)>::value>(
                [&](auto... Is_scal) {
                  // evaluate all expressions
                  this->assignment_impl(std::tuple_cat(
                      nonscalars,
                      this->make_assignment_pair(
                          std::get<2>(std::get<Is_scal>(scalars)),
                          sum_2d(std::get<1>(std::get<Is_scal>(scalars))))...));

                  // copy results from the OpenCL device and increment the
                  // adjoints
                  std::tie(std::get<0>(std::get<Is_scal>(scalars))...)
                      = std::make_tuple(std::get<0>(std::get<Is_scal>(scalars))
                                        + sum(from_matrix_cl(std::get<2>(
                                            std::get<Is_scal>(scalars))))...);
                });
          });
    });
  }

 private:
  /**
   * Selects assignments that have scalar var results.
   * @tparam T_expression type of expression
   * @param result result
   * @param expression expression
   * @return triplet of reference to adjoint, expression and temporary
   * `matrix_cl`
   */
  template <typename T_expression>
  auto select_scalar_assignments(const var& result, T_expression&& expression) {
    return std::make_tuple(std::tuple<double&, T_expression, matrix_cl<double>>(
        result.adj(), std::forward<T_expression>(expression), {}));
  }
  /**
   * Selects assignments that have scalar var results.
   * @tparam T_result type of result. This overload is used for non-scalar-var
   * results.
   * @tparam T_expression type of expression
   * @param result result
   * @param expression expression
   * @return empty tuple
   */
  template <typename T_result, typename T_expression,
            require_all_not_same_t<T_result, var>* = nullptr>
  auto select_scalar_assignments(T_result&& result, T_expression&& expression) {
    return std::make_tuple();
  }

  /**
   * Selects assignments that have non-scalar var results.
   * @tparam T_result type of result. This overload is used for non-scalar vars.
   * @tparam T_expression type of expression
   * @param result result
   * @param expression expression
   * @return pair of result and expression or empty tuple (if the result is
   * check or the expression is `calc_if<false,T>`.
   */
  template <typename T_result, typename T_expression,
            require_not_stan_scalar_t<T_result>* = nullptr,
            require_st_var<T_result>* = nullptr>
  auto select_nonscalar_assignments(const T_result& result,
                                    T_expression&& expression) {
    return results_cl<T_results...>::make_assignment_pair(
        result.adj(), std::forward<T_expression>(expression));
  }
  /**
   * Selects assignments that have non-scalar var results.
   * @tparam T_result type of result. This overload is used for results that are
   * either scalars or not vars.
   * @tparam T_expression type of expression
   * @param result result
   * @param expression expression
   * @return empty tuple
   */
  template <
      typename T_result, typename T_expression,
      std::enable_if_t<is_stan_scalar<T_result>::value
                       || !is_var<scalar_type_t<T_result>>::value>* = nullptr>
  auto select_nonscalar_assignments(T_result&& result,
                                    T_expression&& expression) {
    return std::make_tuple();
  }
};

/**
 * Deduces types for constructing `adjoint_results_cl` object.
 * @tparam T_results types of results
 * @param results results that will be calculated.
 */
template <typename... T_results>
adjoint_results_cl<T_results...> adjoint_results(T_results&&... results) {
  return adjoint_results_cl<T_results...>(std::forward<T_results>(results)...);
}

}  // namespace math
}  // namespace stan

#endif
#endif
