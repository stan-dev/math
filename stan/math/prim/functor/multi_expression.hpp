#ifndef STAN_MATH_PRIM_FUNCTOR_MULTI_EXPRESSION_HPP
#define STAN_MATH_PRIM_FUNCTOR_MULTI_EXPRESSION_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <tuple>

namespace stan {
namespace math {

/**
 * Represents multiple exprs that will be calculated in same kernel.
 * @tparam T_expressions types of exprs
 */
template <typename... T_expressions>
class eigen_expressions_ {
 public:
  /**
   * Constructor.
   * @param exprs exprs that will be calculated in same kernel.
   */
  explicit eigen_expressions_(T_expressions&&... exprs)
      : exprs_(std::forward<T_expressions>(exprs)...) {}

 private:
  std::tuple<T_expressions...> exprs_;
  template <typename... T_results>
  friend class eigen_results_;
};

/**
 * Deduces types for constructing \c expressions_ object.
 * @tparam T_expressions types of exprs
 * @param exprs exprs that will be used in same kernel.
 */
template <typename... T_expressions,
          require_all_eigen_t<T_expressions...>* = nullptr>
eigen_expressions_<T_expressions...> eigen_expressions(
    T_expressions&&... exprs) {
  return eigen_expressions_<T_expressions...>(
      std::forward<T_expressions>(exprs)...);
}

/**
 * Represents results that will be calculated in same loop.
 * @tparam T_results types of results
 */
template <typename... T_results>
class eigen_results_ {
  std::tuple<T_results...> results_;

  /**
   *
   * This functor is called using Eigen's NullaryExpr.
   */
  template <typename... T_expressions>
  class assignment_functor {
    std::tuple<T_results...>& results_;
    const std::tuple<T_expressions...>& exprs_;

   public:
    assignment_functor(std::tuple<T_results...>& results,
                       const std::tuple<T_expressions...>& expressions)
        : results_(results), exprs_(expressions) {}

    inline decltype(auto) operator()(Eigen::Index row, Eigen::Index col) const {
      index_apply<sizeof...(T_results) - 1>([&](auto... Is) {
        static_cast<void>(std::initializer_list<int>{
            (std::get<Is + 1>(results_).coeffRef(row, col)
             = std::get<Is + 1>(exprs_).coeff(row, col),
             0)...});
      });
      return std::get<0>(exprs_).coeff(row, col);
    }

    inline decltype(auto) operator()(Eigen::Index index) const {
      index_apply<sizeof...(T_results) - 1>([&](auto... Is) {
        static_cast<void>(std::initializer_list<int>{
            (std::get<Is + 1>(results_).coeffRef(index)
             = std::get<Is + 1>(exprs_).coeff(index),
             0)...});
      });
      return std::get<0>(exprs_).coeff(index);
    }
  };

 public:
  /**
   * Constructor.
   * @param results results that will be calculated in same kernel
   */
  explicit eigen_results_(T_results&&... results)
      : results_(std::forward<T_results>(results)...) {}

  /**
   * Assigning \c expressions_cl object to \c eigen_results_ object generates
   * and executes the kernel that evaluates expressions and stores them into
   * result expressions this object contains.
   * @tparam T_expressions types of expressions
   * @param expressions expressions
   */
  template <typename... T_expressions,
            typename = std::enable_if_t<sizeof...(T_results)
                                        == sizeof...(T_expressions)>>
  void operator=(const eigen_expressions_<T_expressions...>& expressions) {
    using T_first_expr = std::tuple_element_t<0, std::tuple<T_expressions...>>;
    const auto& first_expr = std::get<0>(expressions.exprs_);
    index_apply<sizeof...(T_results) - 1>([&](auto... Is) {
      static_cast<void>(std::initializer_list<int>{
          (Eigen::internal::resize_if_allowed(
               std::get<Is + 1>(results_), std::get<Is + 1>(expressions.exprs_),
               Eigen::internal::assign_op<int, int>()),
           // types in the assign_op don't matter for the resizing
           0)...});
    });
    std::get<0>(results_) = T_first_expr::NullaryExpr(
        first_expr.rows(), first_expr.cols(),
        assignment_functor<T_expressions...>(results_, expressions.exprs_));
  }
};

/**
 * Deduces types for constructing \c results_cl object.
 * @tparam T_results types of results
 * @param results results that will be calculated in same kernel.
 */
template <typename... T_results, require_all_eigen_t<T_results...>* = nullptr>
eigen_results_<T_results...> eigen_results(T_results&&... results) {
  return eigen_results_<T_results...>(std::forward<T_results>(results)...);
}

}  // namespace math
}  // namespace stan
#endif
