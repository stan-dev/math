#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_consistent_sizes.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/prim/functor/apply.hpp>

#include <cstddef>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace stan {
namespace math {

namespace internal {

template <typename T, typename... Types>
inline void check_all_matching_dims(const char* function, const T& x,
                                    const Types&... xs) {
  std::size_t arg_idx = 2;
  (
      [&](const auto& y) {
        if (x.rows() != y.rows() || x.cols() != y.cols()) {
          [&]() STAN_COLD_PATH {
            const std::string name = "x" + std::to_string(arg_idx);
            check_matching_dims(function, "x1", x, name.c_str(), y);
          }();
        }
        ++arg_idx;
      }(xs),
      ...);
}

}  // namespace internal

/**
 * Apply a functor to each element of a std::vector and collect the results.
 *
 * For a vector `x` of length `n`, returns a std::vector whose `i`-th element
 * is `f(x[i], args...)`. Additional arguments are shared across calls, which
 * allows passing parameters without closures. If `x` is empty, returns an
 * empty std::vector.
 *
 * The element type of the returned std::vector is deduced from the return type
 * of `f` to allow for cases where the input and return scalar types differ
 * (e.g., functions implicitly promoting integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam T Type of input std::vector.
 * @tparam Args Types of shared arguments passed to each call.
 * @param f functor to apply to each element.
 * @param x std::vector input to which operation is applied.
 * @param args shared arguments passed to each call.
 * @return std::vector with result of applying functor to each element of `x`.
 */
template <typename F, typename T, typename... Args,
          require_std_vector_t<T>* = nullptr>
inline auto map(F&& f, T&& x, Args&&... args) {
  using T_return = std::decay_t<decltype(f(x[0], args...))>;
  std::vector<T_return> result;
  result.reserve(x.size());

  for (auto&& xi : x) {
    result.push_back(f(xi, args...));
  }

  return result;
}

/**
 * Apply a functor elementwise to N std::vector inputs of the same size.
 *
 * For vectors `args...` of length `n`, returns a std::vector whose `i`-th
 * element is `f(args[0][i], args[1][i], ...)`. If the inputs are empty,
 * returns an empty std::vector.
 *
 * Unlike `map`, this overload only zips containers and does not accept
 * additional shared trailing arguments. Pass shared state through the
 * functor (e.g. a capturing lambda) when needed.
 *
 * The element type of the returned std::vector is deduced from the return type
 * of `f` to allow for cases where the input and return scalar types differ
 * (e.g., functions implicitly promoting integers).
 *
 * @tparam F Type of functor to apply.
 * @tparam Types Types of input std::vector containers.
 * @param f functor to apply to each tuple of elements.
 * @param args std::vector inputs to which operation is applied.
 * @return std::vector with result of applying functor to each tuple of
 * elements.
 * @throw std::invalid_argument if the inputs have different sizes.
 */
template <typename F, typename... Types,
          require_all_std_vector_t<Types...>* = nullptr>
inline auto mapN(F&& f, Types&&... args) {
  static_assert(sizeof...(Types) >= 2,
                "mapN requires at least two std::vector inputs.");
  static constexpr const char* function = "mapN";
  check_consistent_sizes(function, args...);

  const std::size_t n = std::get<0>(std::forward_as_tuple(args...)).size();
  using T_return = std::decay_t<decltype(f((args[0])...))>;
  std::vector<T_return> result;
  result.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    result.push_back(f((args[i])...));
  }

  return result;
}

/**
 * Apply a functor to each row of an Eigen matrix and collect the results.
 *
 * For a matrix `m` with `n` rows, returns a matrix whose `i`-th row is
 * `f(m.row(i), args...)`. Additional arguments are shared across calls. If
 * `m` has zero rows, returns a 0x0 matrix (output column count is unknown
 * when the functor is never called). If the first returned row is empty,
 * returns a 0x0 matrix without calling the functor on the remaining rows.
 *
 * Expensive Eigen expressions passed as `m` are evaluated once via `to_ref`
 * before the row loop, so compound expressions (e.g. `A * B`) are not
 * re-evaluated for every row.
 *
 * The functor must return a type assignable to a matrix row. All returned
 * rows must have the same number of columns.
 *
 * @tparam F Type of functor to apply.
 * @tparam T Eigen matrix type.
 * @tparam Args Types of shared arguments passed to each call.
 * @param f functor to apply to each row.
 * @param m Eigen matrix input to which operation is applied.
 * @param args shared arguments passed to each call.
 * @return Eigen matrix with result of applying functor to each row of `m`.
 * @throw std::invalid_argument if returned rows have inconsistent sizes.
 */
template <typename F, typename T, typename... Args,
          require_eigen_matrix_dynamic_t<T>* = nullptr>
inline auto row_map(F&& f, T&& m, Args&&... args) {
  static constexpr const char* function = "row_map";
  decltype(auto) m_ref = to_ref(std::forward<T>(m));
  using result_row_t = plain_type_t<decltype(f(m_ref.row(0), args...))>;
  using T_return = scalar_type_t<result_row_t>;
  using matrix_t = Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic>;

  const Eigen::Index n_rows = m_ref.rows();
  if (n_rows == 0) {
    return matrix_t(0, 0);
  }

  result_row_t first_row = f(m_ref.row(0), args...);
  if (first_row.size() == 0) {
    return matrix_t(0, 0);
  }
  matrix_t result(n_rows, first_row.cols());
  result.row(0) = first_row;
  for (Eigen::Index i = 1; i < n_rows; ++i) {
    result_row_t row = f(m_ref.row(i), args...);
    check_size_match(function, "columns of result", result.cols(),
                     "columns of returned row", row.cols());
    result.row(i) = row;
  }
  return result;
}

/**
 * Apply a functor to each column of an Eigen matrix and collect the results.
 *
 * For a matrix `m` with `k` columns, returns a matrix whose `j`-th column is
 * `f(m.col(j), args...)`. Additional arguments are shared across calls. If
 * `m` has zero columns, returns a 0x0 matrix (output row count is unknown
 * when the functor is never called). If the first returned column is empty,
 * returns a 0x0 matrix without calling the functor on the remaining columns.
 *
 * Expensive Eigen expressions passed as `m` are evaluated once via `to_ref`
 * before the column loop, so compound expressions (e.g. `A * B`) are not
 * re-evaluated for every column.
 *
 * The functor must return a type assignable to a matrix column. All returned
 * columns must have the same number of rows.
 *
 * @tparam F Type of functor to apply.
 * @tparam T Eigen matrix type.
 * @tparam Args Types of shared arguments passed to each call.
 * @param f functor to apply to each column.
 * @param m Eigen matrix input to which operation is applied.
 * @param args shared arguments passed to each call.
 * @return Eigen matrix with result of applying functor to each column of `m`.
 * @throw std::invalid_argument if returned columns have inconsistent sizes.
 */
template <typename F, typename T, typename... Args,
          require_eigen_matrix_dynamic_t<T>* = nullptr>
inline auto col_map(F&& f, T&& m, Args&&... args) {
  static constexpr const char* function = "col_map";
  decltype(auto) m_ref = to_ref(std::forward<T>(m));
  using result_col_t = plain_type_t<decltype(f(m_ref.col(0), args...))>;
  using T_return = scalar_type_t<result_col_t>;
  using matrix_t = Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic>;

  const Eigen::Index n_cols = m_ref.cols();
  if (n_cols == 0) {
    return matrix_t(0, 0);
  }

  result_col_t first_col = f(m_ref.col(0), args...);
  if (first_col.size() == 0) {
    return matrix_t(0, 0);
  }
  matrix_t result(first_col.rows(), n_cols);
  result.col(0) = first_col;
  for (Eigen::Index j = 1; j < n_cols; ++j) {
    result_col_t col = f(m_ref.col(j), args...);
    check_size_match(function, "rows of result", result.rows(),
                     "rows of returned column", col.rows());
    result.col(j) = col;
  }
  return result;
}

/**
 * Apply a functor rowwise to N Eigen matrices with identical dimensions.
 *
 * For matrices `args...` with `n` rows, returns a matrix whose `i`-th row is
 * `f(args[0].row(i), args[1].row(i), ...)`. If the inputs have zero rows,
 * returns a 0x0 matrix. If the first returned row is empty, returns a 0x0
 * matrix without calling the functor on the remaining rows.
 *
 * Unlike `row_map`, this overload only zips matrices and does not accept
 * additional shared trailing arguments.
 *
 * Expensive Eigen expressions in `args...` are each evaluated once via
 * `to_ref` before the row loop.
 *
 * The functor must return a type assignable to a matrix row. All returned
 * rows must have the same number of columns.
 *
 * @tparam F Type of functor to apply.
 * @tparam Types Eigen matrix types.
 * @param f functor to apply to each tuple of rows.
 * @param args Eigen matrix inputs to which operation is applied.
 * @return Eigen matrix with result of applying functor to each tuple of rows.
 * @throw std::invalid_argument if the inputs have different dimensions or
 * returned rows have inconsistent sizes.
 */
template <typename F, typename... Types,
          require_all_eigen_matrix_dynamic_t<Types...>* = nullptr>
inline auto row_mapN(F&& f, Types&&... args) {
  static_assert(sizeof...(Types) >= 2,
                "row_mapN requires at least two Eigen matrix inputs.");
  static constexpr const char* function = "row_mapN";
  internal::check_all_matching_dims(function, args...);

  auto m_refs = std::tuple{to_ref(std::forward<Types>(args))...};
  const Eigen::Index n_rows = std::get<0>(m_refs).rows();
  using result_row_t = plain_type_t<decltype(f(
      (std::declval<ref_type_t<Types&&>>().row(0))...))>;
  using T_return = scalar_type_t<result_row_t>;
  using matrix_t = Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic>;

  if (n_rows == 0) {
    return matrix_t(0, 0);
  }

  return apply(
      [&](auto&&... ms) {
        result_row_t first_row = f((ms.row(0))...);
        if (first_row.size() == 0) {
          return matrix_t(0, 0);
        }
        matrix_t result(n_rows, first_row.cols());
        result.row(0) = first_row;
        for (Eigen::Index i = 1; i < n_rows; ++i) {
          result_row_t row = f((ms.row(i))...);
          check_size_match(function, "columns of result", result.cols(),
                           "columns of returned row", row.cols());
          result.row(i) = row;
        }
        return result;
      },
      m_refs);
}

/**
 * Apply a functor columnwise to N Eigen matrices with identical dimensions.
 *
 * For matrices `args...` with `k` columns, returns a matrix whose `j`-th
 * column is `f(args[0].col(j), args[1].col(j), ...)`. If the inputs have
 * zero columns, returns a 0x0 matrix. If the first returned column is empty,
 * returns a 0x0 matrix without calling the functor on the remaining columns.
 *
 * Unlike `col_map`, this overload only zips matrices and does not accept
 * additional shared trailing arguments.
 *
 * Expensive Eigen expressions in `args...` are each evaluated once via
 * `to_ref` before the column loop.
 *
 * The functor must return a type assignable to a matrix column. All returned
 * columns must have the same number of rows.
 *
 * @tparam F Type of functor to apply.
 * @tparam Types Eigen matrix types.
 * @param f functor to apply to each tuple of columns.
 * @param args Eigen matrix inputs to which operation is applied.
 * @return Eigen matrix with result of applying functor to each tuple of
 * columns.
 * @throw std::invalid_argument if the inputs have different dimensions or
 * returned columns have inconsistent sizes.
 */
template <typename F, typename... Types,
          require_all_eigen_matrix_dynamic_t<Types...>* = nullptr>
inline auto col_mapN(F&& f, Types&&... args) {
  static_assert(sizeof...(Types) >= 2,
                "col_mapN requires at least two Eigen matrix inputs.");
  static constexpr const char* function = "col_mapN";
  internal::check_all_matching_dims(function, args...);

  auto m_refs = std::tuple{to_ref(std::forward<Types>(args))...};
  const Eigen::Index n_cols = std::get<0>(m_refs).cols();
  using result_col_t = plain_type_t<decltype(f(
      (std::declval<ref_type_t<Types&&>>().col(0))...))>;
  using T_return = scalar_type_t<result_col_t>;
  using matrix_t = Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic>;
  if (n_cols == 0) {
    return matrix_t(0, 0);
  }
  return apply(
      [&](auto&&... ms) {
        result_col_t first_col = f((ms.col(0))...);
        if (first_col.size() == 0) {
          return matrix_t(0, 0);
        }
        matrix_t result(first_col.rows(), n_cols);
        result.col(0) = first_col;
        for (Eigen::Index j = 1; j < n_cols; ++j) {
          result_col_t col = f((ms.col(j))...);
          check_size_match(function, "rows of result", result.rows(),
                           "rows of returned column", col.rows());
          result.col(j) = col;
        }
        return result;
      },
      m_refs);
}

}  // namespace math
}  // namespace stan

#endif
