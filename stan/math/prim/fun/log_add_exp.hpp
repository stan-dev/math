#ifndef STAN_MATH_PRIM_FUN_LOG_ADD_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG_ADD_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>

namespace stan {
namespace math {

/**
 * Calculates the elementwise sum of exponentials without overflow.
 *
 * \f$\log (\exp(a) + \exp(b)) = m + \log(\exp(a-m) + \exp(b-m))\f$,
 *
 * where \f$m = max(a, b)\f$.
 *
 * @tparam T1 type of the first variable
 * @tparam T2 type of the second variable
 * @param a the first variable
 * @param b the second variable
 */

template <typename T1, typename T2, require_all_not_st_var<T1, T2>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_add_exp(const T2& a, const T1& b) {

  if (a == NEGATIVE_INFTY) {
    return b;
  }
  if (b == NEGATIVE_INFTY) {
    return a;
  }
  if (a == INFTY || b == INFTY) {
    return INFTY;
  }

  const double max_val = std::max(a, b);
  return max_val + std::log(std::exp(a - max_val) + std::exp(b - max_val));
}

/**
 * Calculates the element-wise log sum of exponentials for two containers.
 * For vectors a and b, computes log(exp(a[i]) + exp(b[i])) for each element i.
 * If sizes don't match, uses the smaller size.
 *
 * @tparam T1 type of first container
 * @tparam T2 type of second container
 * @param a First input container
 * @param b Second input container
 * @return Container with element-wise log_add_exp results
 */
template <typename T, require_container_st<std::is_arithmetic, T>* = nullptr>
inline auto log_add_exp(const T& a, const T& b) {

    // Check if sizes are compatible
    if constexpr (stan::is_eigen<T>::value) {
        // Check if both matrices/vectors have the same dimensions
        stan::math::check_matching_dims("log_add_exp", "a", a, "b", b);

        // Determine the number of rows and columns for the result
        size_t rows = a.rows();
        size_t cols = b.cols();
        using return_t = return_type_t<T>;

        Eigen::Matrix<return_t, Eigen::Dynamic, Eigen::Dynamic> result(rows, cols);

        // Iterate over each element
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                double a_val = (a.cols() == 1) ? a(i, 0) : a(i, j); // Handle column vector or matrix
                double b_val = (b.rows() == 1) ? b(0, j) : b(i, j); // Handle row vector or matrix

                if (a_val == NEGATIVE_INFTY) {
                    result(i, j) = b_val;
                } else if (b_val == NEGATIVE_INFTY) {
                    result(i, j) = a_val;
                } else if (a_val == INFTY || b_val == INFTY) {
                    result(i, j) = INFTY;
                } else {
                    result(i, j) = log_sum_exp(a_val, b_val);
                }
            }
        }

        return result;
    } else if constexpr (std::is_same_v<T, std::vector<typename T::value_type>>) {
        // Handle std::vector
        if (a.size() != b.size()) {
            throw std::invalid_argument("Sizes of x and y must match.");
        }

        using return_t = return_type_t<T>;
        std::vector<return_t> result(a.size());

        for (size_t i = 0; i < a.size(); ++i) {
            double a_val = a[i];
            double b_val = b[i];

            if (a_val == NEGATIVE_INFTY) {
                result[i] = b_val;
            } else if (b_val == NEGATIVE_INFTY) {
                result[i] = a_val;
            } else if (a_val == INFTY || b_val == INFTY) {
                result[i] = INFTY;
            } else {
                result[i] = log_sum_exp(a_val, b_val);
            }
        }

        return result;
    } else {
        throw std::invalid_argument("Unsupported container type.");
    }
}

/**
 *  Enables the vectorized application of the log_add_exp function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1
 * @tparam T2
 * @param a
 * @param b
 * @return auto
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_add_exp(const T1& a, const T2& b) {
    // Check if both are Eigen/vectors
    if constexpr (stan::is_eigen<T1>::value && stan::is_eigen<T2>::value) {
        // Check if both matrices/vectors have the same dimensions
        stan::math::check_matching_dims("log_add_exp", "a", a, "b", b);
    } else {
        // Check if sizes are compatible for other types
        if (a.size() != b.size()) {
            throw std::invalid_argument("Sizes of x and y must match or be compatible.");
        }
    }

    // If dimensions are verified to match, apply the operation
    return apply_scalar_binary(
      a, b, [](const auto& c, const auto& d) { return log_add_exp(c, d); });
}

}  // namespace math
}  // nfamespace stan

#endif
