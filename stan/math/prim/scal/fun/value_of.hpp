#ifndef STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_VALUE_OF_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {

  namespace math {

    /**
     * Return the value of the specified scalar argument
     * converted to a double value.
     *
     * <p>See the <code>stan::math::primitive_value</code> function to
     * extract values without casting to <code>double</code>.
     *
     * <p>This function is meant to cover the primitive types. For
     * types requiring pass-by-reference, this template function
     * should be specialized.
     *
     * @tparam T Type of scalar.
     * @param x Scalar to convert to double.
     * @return Value of scalar cast to a double.
     */
    template <typename T>
    inline double value_of(const T x) {
      return static_cast<double>(x);
    }

    /**
     * Return the specified argument.
     *
     * <p>See <code>value_of(T)</code> for a polymorphic
     * implementation using static casts.
     *
     * <p>This inline pass-through no-op should be compiled away.
     *
     * @param x Specified value.
     * @return Specified value.
     */
    template <>
    inline double value_of<double>(const double x) {
      return x;
    }

    template <typename T>
    inline std::vector<double> value_of(const std::vector<T>& x) {
      size_t size = x.size();
      std::vector<double> result(size);
      for (int i=0; i < size; i++)
        result[i] = value_of(x[i]);
      return result;
    }

    template <typename T, int R, int C>
    inline typename Eigen::Matrix<double, R, C>
    value_of(const Eigen::Matrix<T, R, C>& x) {
      int R_ = x.rows();
      int C_ = x.cols();
      int S = x.size();
      Eigen::Matrix<double, R, C> result(R_, C_);
      double* datap = result.data();
      const T* datax = value_of(x.data());
      for (int i=0; i < S; i++)
        datap[i] = value_of(datax[i]);
      return result;
    }

    template <>
    inline std::vector<double> value_of(const std::vector<double>& x) {
      return x;
    }

    template <int R, int C>
    inline typename Eigen::Matrix<double, R, C>
    value_of(const Eigen::Matrix<double, R, C>& x) {
      return x;
    }

  }
}

#endif
