#ifndef STAN_MATH_PRIM_ARR_FUN_APPEND_ARRAY_HPP
#define STAN_MATH_PRIM_ARR_FUN_APPEND_ARRAY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <vector>
#include <iostream>

namespace stan {
  namespace math {
    /**
     * This is custom return type promotion logic for append_array.
     *
     * By default, the type member is computed from the return_type metaprogram.
     *
     * Specializations are provided for std::vector and Eigen::Matrix so that
     * the promotion happens only on the scalar types inside those containers
     * without changing the containers themselves. Roughly,
     * std::vector<T1>, std::vector<T2> promotes to
     * std::vector<return_type<T1, T2>>.
     *
     * @tparam T1 First type to be promoted
     * @tparam T2 Second type to be promoted
     */
    template <typename T1, typename T2>
    struct append_return_type {
      typedef typename return_type<T1, T2>::type type;
    };

    /**
     * Part of return type promotion logic for append_array
     *
     * return_type<int, int> promotes to double, which is not the desired
     * behavior (hence this specialization)
     *
     * @tparam T1 First type to be promoted
     * @tparam T2 Second type to be promoted
     */
    template <>
    struct append_return_type<int, int> {
      typedef int type;
    };

    /**
     * Part of return type promotion logic for append_array
     *
     * The type member is set to an Eigen Matrix with a scalar type that is the
     * result of promoting the scalar types of the two provided Eigen Matrix
     * types
     *
     * @tparam T1 Scalar type of first matrix type
     * @tparam T2 Scalar type of first matrix type
     * @tparam R Eigen RowsAtCompileTime of both matrices
     * @tparam C Eigen ColsAtCompileTime of both matrices
     */
    template <typename T1, typename T2, int R, int C>
    struct append_return_type<Eigen::Matrix<T1, R, C>,
                              Eigen::Matrix<T2, R, C> > {
      typedef typename
      Eigen::Matrix<typename append_return_type<T1, T2>::type, R, C> type;
    };

    /**
     * Part of the return type promotion logic for append_array
     *
     * If the types of both template arguments are std::vectors, the type member
     * is recursively computed as the append_return_type of the scalar types
     * associated with those std::vectors.
     */
    template <typename T1, typename T2>
    struct append_return_type<std::vector<T1>, std::vector<T2> > {
      typedef typename
      std::vector<typename append_return_type<T1, T2>::type> type;
    };

    /**
     * Part of copy logic for append_array.
     *
     * This copies the value at x to y. By default, it is assumed the necessary
     * copy constructors are available to do the casting from T1 to T2.
     *
     * Specializations for std::vector and Eigen::Matrix are provided to do
     * copy with casting according to the types computed with append_return_type
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
     * Part of copy logic for append_array.
     *
     * The necessary copying/casting for Eigen matrices is delegated to Eigen
     * template cast
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
     * Part of copy logic for append_array.
     *
     * This specialization deals with copying vectors. First, y is resized to
     * match the size of x, and then append_array_copy is recursively called on
     * each element of x to copy it to an element of y.
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
     * Part of size checking logic for append_array
     *
     * The goal of these functions is to check if the two arguments have
     * compatible sizes for append_array_copy. The default behavior is to
     * assume that sizes are compatible. Specializations are defined for cases
     * where this is false.
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
     * Part of size checking logic for append_array
     *
     * This specialization deals with Eigen::Matrices, and delegates to
     * check_matching_dims to make sure the two arguments are of compatible
     * sizes.
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
     * Part of custom size checking logic for append_array
     *
     * This specialization deals with std::vectors. This function checks that
     * the size of x and y are the same and then recursively checks that the
     * elements are of compatible sizes as well.
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
