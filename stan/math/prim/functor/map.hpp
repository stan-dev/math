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
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/**
 * Return a std::vector whose i-th element is `f(x[i], args...)`.
 * Additional arguments are shared across calls.  If `x` is empty, the result
 * is empty.  The element type of the result is deduced from `f`.
 *
 * @tparam F type of functor
 * @tparam T type of input std::vector
 * @tparam Args types of shared arguments
 * @param[in] f functor applied to each element
 * @param[in] x input std::vector
 * @param[in] args shared arguments passed to each call
 * @return std::vector of results
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
 * Return a std::vector whose i-th element is
 * `f(get<0>(xs)[i], get<1>(xs)[i], ..., args...)`.
 * The tuple `xs` must contain at least two equal-sized std::vector inputs.
 * Additional trailing arguments are shared across calls.  If the inputs are
 * empty, the result is empty.  The element type of the result is deduced from
 * `f`.
 *
 * @tparam F type of functor
 * @tparam Tuple type of tuple of std::vector inputs
 * @tparam Args types of shared arguments
 * @param[in] f functor applied to each tuple of elements
 * @param[in] xs tuple of std::vector inputs
 * @param[in] args shared arguments passed to each call
 * @return std::vector of results
 * @throw std::invalid_argument if the inputs have different sizes
 */
template <typename F, typename Tuple, typename... Args,
          require_tuple_t<Tuple>* = nullptr>
inline auto mapN(F&& f, Tuple&& xs, Args&&... args) {
  static_assert(tuple_size_v<Tuple> >= 2,
                "mapN requires a tuple of at least two std::vector inputs.");
  static constexpr const char* function = "mapN";

  return apply(
      [&](auto&&... xs) {
        static_assert((is_std_vector<std::decay_t<decltype(xs)>>::value && ...),
                      "mapN tuple elements must be std::vector.");
        check_consistent_sizes(function, xs...);

        const std::size_t n = std::get<0>(std::forward_as_tuple(xs...)).size();
        using T_return = std::decay_t<decltype(f((xs[0])..., args...))>;
        std::vector<T_return> result;
        result.reserve(n);

        for (std::size_t i = 0; i < n; ++i) {
          result.push_back(f((xs[i])..., args...));
        }
        return result;
      },
      std::forward<Tuple>(xs));
}

/**
 * Return a std::vector whose i-th element is `f(x1[i], x2[i], ...)`.
 * Convenience overload equivalent to
 * `mapN(f, std::forward_as_tuple(x1, x2, xs...))`.
 *
 * @tparam F type of functor
 * @tparam T1 type of first std::vector
 * @tparam T2 type of second std::vector
 * @tparam Ts types of additional std::vector inputs
 * @param[in] f functor applied to each tuple of elements
 * @param[in] x1 first std::vector
 * @param[in] x2 second std::vector
 * @param[in] xs additional std::vector inputs
 * @return std::vector of results
 * @throw std::invalid_argument if the inputs have different sizes
 */
template <typename F, typename T1, typename T2, typename... Ts,
          require_all_std_vector_t<T1, T2, Ts...>* = nullptr>
inline auto mapN(F&& f, T1&& x1, T2&& x2, Ts&&... xs) {
  return mapN(std::forward<F>(f), std::forward_as_tuple(x1, x2, xs...));
}

/**
 * Return a matrix whose i-th row is `f(m.row(i), args...)`.
 * Additional arguments are shared across calls.  If `m` has zero rows, or if
 * the first returned row is empty, returns a 0x0 matrix.  All returned rows
 * must have the same number of columns.
 *
 * @tparam F type of functor
 * @tparam T type of Eigen matrix
 * @tparam Args types of shared arguments
 * @param[in] f functor applied to each row
 * @param[in] m input matrix
 * @param[in] args shared arguments passed to each call
 * @return matrix of results
 * @throw std::invalid_argument if returned rows have inconsistent sizes
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
 * Return a matrix whose j-th column is `f(m.col(j), args...)`.
 * Additional arguments are shared across calls.  If `m` has zero columns, or
 * if the first returned column is empty, returns a 0x0 matrix.  All returned
 * columns must have the same number of rows.
 *
 * @tparam F type of functor
 * @tparam T type of Eigen matrix
 * @tparam Args types of shared arguments
 * @param[in] f functor applied to each column
 * @param[in] m input matrix
 * @param[in] args shared arguments passed to each call
 * @return matrix of results
 * @throw std::invalid_argument if returned columns have inconsistent sizes
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
 * Return a matrix whose i-th row is
 * `f(get<0>(ms).row(i), get<1>(ms).row(i), ..., args...)`.
 * The tuple `ms` must contain at least two Eigen matrices with identical
 * dimensions.  Additional trailing arguments are shared across calls.  If the
 * inputs have zero rows, or if the first returned row is empty, returns a 0x0
 * matrix.  All returned rows must have the same number of columns.
 *
 * @tparam F type of functor
 * @tparam Tuple type of tuple of Eigen matrix inputs
 * @tparam Args types of shared arguments
 * @param[in] f functor applied to each tuple of rows
 * @param[in] ms tuple of Eigen matrix inputs
 * @param[in] args shared arguments passed to each call
 * @return matrix of results
 * @throw std::invalid_argument if the inputs have different dimensions or
 * returned rows have inconsistent sizes
 */
template <typename F, typename Tuple, typename... Args,
          require_tuple_t<Tuple>* = nullptr>
inline auto row_mapN(F&& f, Tuple&& ms, Args&&... args) {
  static_assert(tuple_size_v<Tuple> >= 2,
                "row_mapN requires a tuple of at least two Eigen matrix "
                "inputs.");
  static constexpr const char* function = "row_mapN";

  return apply(
      [&](auto&&... mats) {
        static_assert(
            (is_eigen_matrix_dynamic<std::decay_t<decltype(mats)>>::value
             && ...),
            "row_mapN tuple elements must be Eigen matrices.");
        check_matching_dims(function, mats...);

        auto m_refs = std::tuple{to_ref(std::forward<decltype(mats)>(mats))...};
        const Eigen::Index n_rows = std::get<0>(m_refs).rows();
        using result_row_t = plain_type_t<decltype(f(
            (std::declval<ref_type_t<decltype(mats)>>().row(0))..., args...))>;
        using T_return = scalar_type_t<result_row_t>;
        using matrix_t
            = Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic>;

        if (n_rows == 0) {
          return matrix_t(0, 0);
        }

        return apply(
            [&](auto&&... mrs) {
              result_row_t first_row = f((mrs.row(0))..., args...);
              if (first_row.size() == 0) {
                return matrix_t(0, 0);
              }
              matrix_t result(n_rows, first_row.cols());
              result.row(0) = first_row;
              for (Eigen::Index i = 1; i < n_rows; ++i) {
                result_row_t row = f((mrs.row(i))..., args...);
                check_size_match(function, "columns of result", result.cols(),
                                 "columns of returned row", row.cols());
                result.row(i) = row;
              }
              return result;
            },
            m_refs);
      },
      std::forward<Tuple>(ms));
}

/**
 * Return a matrix whose i-th row is `f(m1.row(i), m2.row(i), ...)`.
 * Convenience overload equivalent to
 * `row_mapN(f, std::forward_as_tuple(m1, m2, ms...))`.
 *
 * @tparam F type of functor
 * @tparam T1 type of first Eigen matrix
 * @tparam T2 type of second Eigen matrix
 * @tparam Ts types of additional Eigen matrix inputs
 * @param[in] f functor applied to each tuple of rows
 * @param[in] m1 first Eigen matrix
 * @param[in] m2 second Eigen matrix
 * @param[in] ms additional Eigen matrix inputs
 * @return matrix of results
 * @throw std::invalid_argument if the inputs have different dimensions or
 * returned rows have inconsistent sizes
 */
template <typename F, typename T1, typename T2, typename... Ts,
          require_all_eigen_matrix_dynamic_t<T1, T2, Ts...>* = nullptr>
inline auto row_mapN(F&& f, T1&& m1, T2&& m2, Ts&&... ms) {
  return row_mapN(std::forward<F>(f), std::forward_as_tuple(m1, m2, ms...));
}

/**
 * Return a matrix whose j-th column is
 * `f(get<0>(ms).col(j), get<1>(ms).col(j), ..., args...)`.
 * The tuple `ms` must contain at least two Eigen matrices with identical
 * dimensions.  Additional trailing arguments are shared across calls.  If the
 * inputs have zero columns, or if the first returned column is empty, returns
 * a 0x0 matrix.  All returned columns must have the same number of rows.
 *
 * @tparam F type of functor
 * @tparam Tuple type of tuple of Eigen matrix inputs
 * @tparam Args types of shared arguments
 * @param[in] f functor applied to each tuple of columns
 * @param[in] ms tuple of Eigen matrix inputs
 * @param[in] args shared arguments passed to each call
 * @return matrix of results
 * @throw std::invalid_argument if the inputs have different dimensions or
 * returned columns have inconsistent sizes
 */
template <typename F, typename Tuple, typename... Args,
          require_tuple_t<Tuple>* = nullptr>
inline auto col_mapN(F&& f, Tuple&& ms, Args&&... args) {
  static_assert(tuple_size_v<Tuple> >= 2,
                "col_mapN requires a tuple of at least two Eigen matrix "
                "inputs.");
  static constexpr const char* function = "col_mapN";

  return apply(
      [&](auto&&... mats) {
        static_assert(
            (is_eigen_matrix_dynamic<std::decay_t<decltype(mats)>>::value
             && ...),
            "col_mapN tuple elements must be Eigen matrices.");
        check_matching_dims(function, mats...);

        auto m_refs = std::tuple{to_ref(std::forward<decltype(mats)>(mats))...};
        const Eigen::Index n_cols = std::get<0>(m_refs).cols();
        using result_col_t = plain_type_t<decltype(f(
            (std::declval<ref_type_t<decltype(mats)>>().col(0))..., args...))>;
        using T_return = scalar_type_t<result_col_t>;
        using matrix_t
            = Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic>;

        if (n_cols == 0) {
          return matrix_t(0, 0);
        }

        return apply(
            [&](auto&&... mrs) {
              result_col_t first_col = f((mrs.col(0))..., args...);
              if (first_col.size() == 0) {
                return matrix_t(0, 0);
              }
              matrix_t result(first_col.rows(), n_cols);
              result.col(0) = first_col;
              for (Eigen::Index j = 1; j < n_cols; ++j) {
                result_col_t col = f((mrs.col(j))..., args...);
                check_size_match(function, "rows of result", result.rows(),
                                 "rows of returned column", col.rows());
                result.col(j) = col;
              }
              return result;
            },
            m_refs);
      },
      std::forward<Tuple>(ms));
}

/**
 * Return a matrix whose j-th column is `f(m1.col(j), m2.col(j), ...)`.
 * Convenience overload equivalent to
 * `col_mapN(f, std::forward_as_tuple(m1, m2, ms...))`.
 *
 * @tparam F type of functor
 * @tparam T1 type of first Eigen matrix
 * @tparam T2 type of second Eigen matrix
 * @tparam Ts types of additional Eigen matrix inputs
 * @param[in] f functor applied to each tuple of columns
 * @param[in] m1 first Eigen matrix
 * @param[in] m2 second Eigen matrix
 * @param[in] ms additional Eigen matrix inputs
 * @return matrix of results
 * @throw std::invalid_argument if the inputs have different dimensions or
 * returned columns have inconsistent sizes
 */
template <typename F, typename T1, typename T2, typename... Ts,
          require_all_eigen_matrix_dynamic_t<T1, T2, Ts...>* = nullptr>
inline auto col_mapN(F&& f, T1&& m1, T2&& m2, Ts&&... ms) {
  return col_mapN(std::forward<F>(f), std::forward_as_tuple(m1, m2, ms...));
}

}  // namespace math
}  // namespace stan

#endif
