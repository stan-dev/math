#ifndef STAN_MATH_PRIM_MAT_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_MAT_FUN_VALUE_OF_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/child_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {

    /**
     * Convert a matrix of type T to a matrix of doubles.
     *
     * T must implement value_of. See
     * test/math/fwd/mat/fun/value_of.cpp for fvar and var usage.
     *
     * @tparam T Scalar type in matrix
     * @tparam R Rows of matrix
     * @tparam C Columns of matrix
     * @param[in] M Matrix to be converted
     * @return Matrix of values
     **/
    template <typename T, int R, int C>
    inline typename Eigen::Matrix<typename child_type<T>::type, R, C>
    value_of(const Eigen::Matrix<T, R, C>& M) {
      int R_ = M.rows();
      int C_ = M.cols();
      int S = M.size();
      Eigen::Matrix<double, R, C> result(R_, C_);
      typename child_type<T>::type * datar = result.data();
      const T * dataM = M.data();
      for (int i=0; i < S; i++)
        datar[i] = value_of(dataM[i]);
      return result;
    }

    /**
     * Convert a std::vector of type T to a std::vector of doubles.
     *
     * T must implement value_of. See
     * test/math/fwd/mat/fun/value_of.cpp for fvar and var usage.
     *
     * @tparam T Scalar type in std::vector
     * @param[in] x std::vector to be converted
     * @return std::vector of values
     **/
    template <typename T>
    inline std::vector<typename child_type<T>::type>
    value_of(const std::vector<T>& x) {
      size_t size = x.size();
      std::vector<typename child_type<T>::type> result(size);
      for (int i=0; i < size; i++)
        result[i] = value_of(x[i]);
      return result;
    }

    /**
     * Return the specified argument.
     *
     * <p>See <code>value_of(T)</code> for a polymorphic
     * implementation using static casts.
     *
     * <p>This inline pass-through no-op should be compiled away.
     *
     * @param x Specified std::vector.
     * @return Specified std::vector.
     */
    template <>
    inline std::vector<double> value_of(const std::vector<double>& x) {
      return x;
    }

    /**
     * Return the specified argument.
     *
     * <p>See <code>value_of(T)</code> for a polymorphic
     * implementation using static casts.
     *
     * <p>This inline pass-through no-op should be compiled away.
     *
     * @param x Specified matrix.
     * @return Specified matrix.
     */
    template <int R, int C>
    inline typename Eigen::Matrix<double, R, C>
    value_of(const Eigen::Matrix<double, R, C>& x) {
      return x;
    }

  }
}

#endif
