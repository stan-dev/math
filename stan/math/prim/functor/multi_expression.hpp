#ifndef STAN_MATH_PRIM_FUNCTOR_MULTI_EXPRESSION_HPP
#define STAN_MATH_PRIM_FUNCTOR_MULTI_EXPRESSION_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>
#include <tuple>

namespace stan {
namespace math {

namespace internal {
constexpr int constexpr_sum() { return 0; }

template <typename Arg0, typename... Args>
constexpr int constexpr_sum(Arg0 arg0, Args... args) {
  return arg0 + constexpr_sum(args...);
}

}  // namespace internal

/**
 * Represents multiple exprs that will be calculated in same loop.
 * @tparam T_expressions types of exprs
 */
template <typename... T_expressions>
class eigen_expressions_ {
 public:
  /**
   * Constructor.
   * @param exprs expresions that will be calculated in same loop.
   * @throw invalid_argument expressions have different sizes (for row major
   * expressions rows and columns are swapped for the check)
   */
  explicit eigen_expressions_(T_expressions&&... exprs)
      : exprs_(std::forward<T_expressions>(exprs)...) {
    index_apply<sizeof...(T_expressions) - 1>([&](auto... Is) {
      constexpr auto first_flags = Eigen::internal::evaluator<std::decay_t<
          std::tuple_element_t<0, std::tuple<T_expressions...>>>>::Flags;
      static_cast<void>(std::initializer_list<int>{
          ((((Eigen::internal::evaluator<std::decay_t<std::tuple_element_t<
                  Is + 1, std::tuple<T_expressions...>>>>::Flags
              ^ first_flags)
             & Eigen::RowMajorBit)
                ? check_matching_dims("eigen_expressions_.eigen_expressions_",
                                      "first expression", std::get<0>(exprs_),
                                      "transposed expression",
                                      std::get<Is + 1>(exprs_).transpose())
                : check_matching_dims("eigen_expressions_.eigen_expressions_",
                                      "first expression", std::get<0>(exprs_),
                                      "expression", std::get<Is + 1>(exprs_))),
           0)...});
    });
  }

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
   * Assign expressions to results using linear indexing.
   * @tparam Linear whether to use linear indexing
   * @tparam T_res_evals types of result evaluators
   * @tparam T_expr_evals types of expression evaluators
   * @param res_evals evaluators for results to assign to
   * @param expr_evals evaluators for expressions to assign
   * @param rows number of rows
   * @param cols number of cols
   */
  template <bool Linear, typename... T_res_evals, typename... T_expr_evals,
            std::enable_if_t<Linear>* = nullptr>
  inline void assign(std::tuple<T_res_evals...>& res_evals,
                     const std::tuple<T_expr_evals...>& expr_evals,
                     Eigen::Index rows, Eigen::Index cols) {
    for (size_t i = 0; i < rows * cols; i++) {
      index_apply<sizeof...(T_results)>([&](auto... Is) {
        static_cast<void>(
            std::initializer_list<int>{(std::get<Is>(res_evals).coeffRef(i)
                                        = std::get<Is>(expr_evals).coeff(i),
                                        0)...});
      });
    }
  }
  /**
   * Assign expressions to results using 2d indexing.
   * @tparam Linear whether to use linear indexing
   * @tparam T_res_evals types of result evaluators
   * @tparam T_expr_evals types of expression evaluators
   * @param res_evals evaluators for results to assign to
   * @param expr_evals evaluators for expressions to assign
   * @param rows number of rows
   * @param cols number of cols
   */
  template <bool Linear, typename... T_res_evals, typename... T_expr_evals,
            std::enable_if_t<!Linear>* = nullptr>
  inline void assign(std::tuple<T_res_evals...>& res_evals,
                     const std::tuple<T_expr_evals...>& expr_evals,
                     Eigen::Index rows, Eigen::Index cols) {
    constexpr bool is_first_row_major
        = std::decay_t<decltype(std::get<0>(expr_evals))>::Flags
          & Eigen::RowMajorBit;
    const Eigen::Index outer_dimension = is_first_row_major ? rows : cols;
    const Eigen::Index inner_dimension = is_first_row_major ? cols : rows;
    for (size_t i = 0; i < outer_dimension; i++) {
      for (size_t j = 0; j < inner_dimension; j++) {
        index_apply<sizeof...(T_results)>([&](auto... Is) {
          static_cast<void>(std::initializer_list<int>{
              ((std::decay_t<decltype(std::get<0>(expr_evals))>::Flags
                        & Eigen::RowMajorBit
                    ? std::get<Is>(res_evals).coeffRef(i, j)
                      = std::get<Is>(expr_evals).coeff(i, j)
                    : std::get<Is>(res_evals).coeffRef(j, i)
                      = std::get<Is>(expr_evals).coeff(j, i)),
               0)...});
        });
      }
    }
  }

  /**
   * Selects and calls appropriate `assign`.
   * @tparam T_expressions types of expressions
   * @param expressions expressions
   */
  template <typename... T_expressions,
            std::enable_if_t<sizeof...(T_results)
                             == sizeof...(T_expressions)>* = nullptr>
  void assign_select(const eigen_expressions_<T_expressions...>& expressions) {
    constexpr bool all_linear = std::min(
        {static_cast<bool>(
             Eigen::internal::evaluator<std::decay_t<T_expressions>>::Flags
             & Eigen::LinearAccessBit)...,
         static_cast<bool>(
             Eigen::internal::evaluator<std::decay_t<T_results>>::Flags
             & Eigen::LinearAccessBit)...});
    constexpr int N_row_major = internal::constexpr_sum(
        static_cast<bool>(
            Eigen::internal::evaluator<std::decay_t<T_expressions>>::Flags
            & Eigen::RowMajorBit)...,
        static_cast<bool>(
            Eigen::internal::evaluator<std::decay_t<T_results>>::Flags
            & Eigen::RowMajorBit)...);
    constexpr int N_col_major
        = sizeof...(T_results) + sizeof...(T_expressions) - N_row_major;

    index_apply<sizeof...(T_results)>([&](auto... Is) {
      std::tuple<Eigen::internal::evaluator<
          std::decay_t<decltype(std::get<Is>(results_))>>...>
      result_evaluators(std::get<Is>(results_)...);
      std::tuple<Eigen::internal::evaluator<
          std::decay_t<decltype(std::get<Is>(expressions.exprs_))>>...>
      expression_evaluators(std::get<Is>(expressions.exprs_)...);

      assign<all_linear && (N_col_major == 0 || N_row_major == 0)>(
          result_evaluators, expression_evaluators,
          std::get<0>(expressions.exprs_).rows(),
          std::get<0>(expressions.exprs_).cols());
    });
  }

 public:
  /**
   * Constructor.
   * @param results results that will be calculated in same kernel
   */
  explicit eigen_results_(T_results&&... results)
      : results_(std::forward<T_results>(results)...) {}

  /**
   * Assigning \c expressions_cl object to \c eigen_results_ object evals the
   * expressions into results.
   * @tparam T_expressions types of expressions
   * @param expressions expressions
   */
  template <typename... T_expressions,
            std::enable_if_t<sizeof...(T_results)
                             == sizeof...(T_expressions)>* = nullptr>
  void operator=(const eigen_expressions_<T_expressions...>& expressions) {
    index_apply<sizeof...(T_results)>([&](auto... Is) {
      static_cast<void>(std::initializer_list<int>{
          (Eigen::internal::resize_if_allowed(
               std::get<Is>(results_), std::get<Is>(expressions.exprs_),
               Eigen::internal::assign_op<int, int>()),
           // types in the assign_op don't matter for the resizing
           0)...});
    });

    assign_select(expressions);
  }

  /**
   * Add \c eigen_results_ to \c expressions_cl in place.
   * @tparam T_expressions types of expressions
   * @param expressions expressions
   */
  template <typename... T_expressions,
            std::enable_if_t<sizeof...(T_results)
                             == sizeof...(T_expressions)>* = nullptr>
  void operator+=(const eigen_expressions_<T_expressions...>& expressions) {
    index_apply<sizeof...(T_results)>([&](auto... Is) {
      assign_select(eigen_expressions(
          (std::get<Is>(results_) + std::get<Is>(expressions.exprs_))...));
    });
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
