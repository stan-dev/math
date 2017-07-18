#ifndef STAN_MATH_PRIM_MAT_FUN_APPEND_ARRAY_HPP
#define STAN_MATH_PRIM_MAT_FUN_APPEND_ARRAY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/append_return_type.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <vector>

namespace stan {
  namespace math {
    /**
     * This function copies the value at x to y, casting from the type of T1 to
     * type T2. This implementation assumes the necessary copy constructors are
     * available to do the casting from T1 to T2.
     *
     * Specializations for std::vector and Eigen::Matrix are provided to copy
     * and cast the scalar elements.
     *
     * @tparam Type of source
     * @tparam Type of destination
     * @param x Source variable
     * @param y Target variable
     */
    template<typename T1, typename T2>
    inline void append_array_copy(const T1& x, T2& y) {
      y = T2(x);
    }

    /**
     * This copies and casts the elements of x to y (by delegating to the Eigen
     * template cast)
     *
     * @tparam Scalar type of source matrix
     * @tparam Scalar Type of destination matrix
     * @tparam R Eigen RowsAtCompileTime of both matrices
     * @tparam C Eigen ColsAtCompileTime of both matrices
     * @param x Source matrix
     * @param y Target matrix
     */
    template<typename T1, typename T2, int R, int C>
    inline void append_array_copy(const Eigen::Matrix<T1, R, C>& x,
                     Eigen::Matrix<T2, R, C>& y) {
      y = x.template cast<T2>();
    }

    /**
     * This function copies and casts the elements of std::vector x to std::vector y.
     * First, y is resized to match the size of x, and then append_array_copy
     * is recursively called on each element of x to copy it to an element of y.
     *
     * @tparam Type of source
     * @tparam Type of destination
     * @param x Source variable
     * @param y Target variable
     */
    template<typename T1, typename T2>
    inline void append_array_copy(const std::vector<T1> &x,
                                  std::vector<T2> &y) {
      y.resize(x.size());
      for (size_t i = 0; i < x.size(); i++) {
        append_array_copy(x[i], y[i]);
      }
    }

    /**
     * This function checks the sizes of both arguments are consistent.
     *
     * The basic version assumes variables are of compatible sizes. There are
     * specializations for Eigen::Matrices and std::vectors (where this
     * assumption is not true).
     *
     * @tparam T1 Type of x
     * @tparam T2 Type of y
     * @param x First argument to check
     * @param y Second argument to check
     */
    template<typename T1, typename T2>
    inline void append_array_check_size(const T1& x, const T2& y) {
    }

    /**
     * This function checks if two Eigen::Matrices are of compatible sizes
     * by delegating to check_matching_dims.
     *
     * @tparam T1 Scalar type of x
     * @tparam T2 Scalar type of y
     * @tparam R Eigen RowsAtCompileTime of both matrices
     * @tparam C Eigen ColsAtCompileTime of both matrices
     * @param x First matrix
     * @param y Second matrix
     */
    template<typename T1, typename T2, int R, int C>
    inline void append_array_check_size(const Eigen::Matrix<T1, R, C>& x,
                     const Eigen::Matrix<T2, R, C>& y) {
      check_matching_dims("append_array", "x", x, "y", y);
    }

    /**
     * This function checks if two std::vectors are of compatible sizes and
     * then recursively calls append_array_check_size to see if the first
     * elements of each std::vector have compatible sizes.
     *
     * @tparam T1 Element type of x
     * @tparam T2 Element type of y
     * @param x First std::vector
     * @param y Second std::vector
     */
    template<typename T1, typename T2>
    inline void append_array_check_size(const std::vector<T1> &x,
                                        const std::vector<T2> &y) {
      check_consistent_sizes("append_array", "size of x", x, "size of y", y);
      if (!x.empty() && !y.empty()) {
        append_array_check_size(x[0], y[0]);
      }
    }

    /**
     * Return the concatenation of two specified vectors in the order of
     *   the arguments.
     *
     * @tparam T1 Type of first vector
     * @tparam T2 Type of second vector
     * @param x First vector
     * @param y Second vector
     * @return A vector of x and y concatenated together (in that order)
     */
    template <typename T1, typename T2>
    inline typename append_return_type<T1, T2>::type
    append_array(const T1& x, const T2& y) {
      typename append_return_type<T1, T2>::type z;
      if (!x.empty() && !y.empty()) {
        append_array_check_size(x[0], y[0]);
      }
      z.resize(x.size() + y.size());
      for (size_t i = 0; i < x.size(); i++)
        append_array_copy(x[i], z[i]);
      for (size_t i = 0; i < y.size(); i++)
        append_array_copy(y[i], z[x.size() + i]);
      return z;
    }
  }
}
#endif
