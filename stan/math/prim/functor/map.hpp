#ifndef STAN_MATH_PRIM_FUNCTOR_MAP_HPP
#define STAN_MATH_PRIM_FUNCTOR_MAP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/eval.hpp>
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
          require_all_tuple_elements_t<is_std_vector, Tuple>* = nullptr>
inline auto mapN(F&& f, Tuple&& xs, Args&&... args) {
  static constexpr const char* function = "mapN";

  return apply(
      [&f, &args...](auto&&... xs_inner) {
        check_consistent_sizes(function, xs_inner...);

        const std::size_t n
            = std::get<0>(std::forward_as_tuple(xs_inner...)).size();
        using T_return = std::decay_t<decltype(f((xs_inner[0])..., args...))>;
        std::vector<T_return> result;
        result.reserve(n);

        for (std::size_t i = 0; i < n; ++i) {
          result.push_back(f((xs_inner[i])..., args...));
        }
        return result;
      },
      std::forward<Tuple>(xs));
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

  using matrix_t = Eigen::Matrix<
      scalar_type_t<plain_type_t<decltype(f(m_ref.row(0), args...))>>,
      Eigen::Dynamic, Eigen::Dynamic>;

  const Eigen::Index n_rows = m_ref.rows();
  if (n_rows == 0) {
    return matrix_t(0, 0);
  }

  auto row_i = eval(f(m_ref.row(0), args...));

  if (row_i.size() == 0) {
    return matrix_t(0, 0);
  }
  matrix_t result(n_rows, row_i.cols());
  result.row(0) = row_i;
  for (Eigen::Index i = 1; i < n_rows; ++i) {
    row_i = f(m_ref.row(i), args...);
    check_size_match(function, "columns of result", result.cols(),
                     "columns of returned row", row_i.cols());
    result.row(i) = row_i;
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

  using matrix_t = Eigen::Matrix<
      scalar_type_t<plain_type_t<decltype(f(m_ref.col(0), args...))>>,
      Eigen::Dynamic, Eigen::Dynamic>;

  const Eigen::Index n_cols = m_ref.cols();
  if (n_cols == 0) {
    return matrix_t(0, 0);
  }

  auto col_j = eval(f(m_ref.col(0), args...));

  if (col_j.size() == 0) {
    return matrix_t(0, 0);
  }
  matrix_t result(col_j.rows(), n_cols);
  result.col(0) = col_j;
  for (Eigen::Index j = 1; j < n_cols; ++j) {
    col_j = f(m_ref.col(j), args...);
    check_size_match(function, "rows of result", result.rows(),
                     "rows of returned column", col_j.rows());
    result.col(j) = col_j;
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
  static constexpr const char* function = "row_mapN";

  return apply(
      [&f, &args...](auto&&... mats) {
        static_assert(
            (is_eigen_matrix_dynamic<std::decay_t<decltype(mats)>>::value
             && ...),
            "row_mapN tuple elements must be Eigen matrices.");
        check_matching_dims(function, mats...);

        auto m_refs = std::tuple{to_ref(std::forward<decltype(mats)>(mats))...};
        const Eigen::Index n_rows = std::get<0>(m_refs).rows();

        return apply(
            [&f, &args..., &n_rows](auto&&... mrs) {
              using matrix_t
                  = Eigen::Matrix<scalar_type_t<plain_type_t<decltype(f(
                                      (mrs.row(0))..., args...))>>,
                                  Eigen::Dynamic, Eigen::Dynamic>;

              if (n_rows == 0) {
                return matrix_t(0, 0);
              }

              auto row_i = eval(f((mrs.row(0))..., args...));

              if (row_i.size() == 0) {
                return matrix_t(0, 0);
              }
              matrix_t result(n_rows, row_i.cols());
              result.row(0) = row_i;
              for (Eigen::Index i = 1; i < n_rows; ++i) {
                row_i = f((mrs.row(i))..., args...);
                check_size_match(function, "columns of result", result.cols(),
                                 "columns of returned row", row_i.cols());
                result.row(i) = row_i;
              }
              return result;
            },
            m_refs);
      },
      std::forward<Tuple>(ms));
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
  static constexpr const char* function = "col_mapN";

  return apply(
      [&f, &args...](auto&&... mats) {
        static_assert(
            (is_eigen_matrix_dynamic<std::decay_t<decltype(mats)>>::value
             && ...),
            "col_mapN tuple elements must be Eigen matrices.");
        check_matching_dims(function, mats...);

        auto m_refs = std::tuple{to_ref(std::forward<decltype(mats)>(mats))...};
        const Eigen::Index n_cols = std::get<0>(m_refs).cols();

        return apply(
            [&f, &args..., &n_cols](auto&&... mrs) {
              using matrix_t
                  = Eigen::Matrix<scalar_type_t<plain_type_t<decltype(f(
                                      (mrs.col(0))..., args...))>>,
                                  Eigen::Dynamic, Eigen::Dynamic>;

              if (n_cols == 0) {
                return matrix_t(0, 0);
              }

              auto col_j = eval(f((mrs.col(0))..., args...));

              if (col_j.size() == 0) {
                return matrix_t(0, 0);
              }
              matrix_t result(col_j.rows(), n_cols);
              result.col(0) = col_j;
              for (Eigen::Index j = 1; j < n_cols; ++j) {
                col_j = f((mrs.col(j))..., args...);
                check_size_match(function, "rows of result", result.rows(),
                                 "rows of returned column", col_j.rows());
                result.col(j) = col_j;
              }
              return result;
            },
            m_refs);
      },
      std::forward<Tuple>(ms));
}

}  // namespace math
}  // namespace stan

#endif
