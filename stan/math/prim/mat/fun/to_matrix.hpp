#ifndef STAN_MATH_PRIM_MAT_FUN_TO_MATRIX_HPP
#define STAN_MATH_PRIM_MAT_FUN_TO_MATRIX_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

namespace stan {
  namespace math {
    /**
     * Returns a matrix with dynamic dimensions constructed from the data passed
     * in. This could be a row vector, column vector, or matrix. The runtime
     * dimensions will be the same as the input.
     *
     * @tparam T type of the scalar
     * @tparam R number of rows
     * @tparam C number of columns
     * @param matrix matrix
     * @return the matrix representation of the input
     */
    template <typename T, int R, int C>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(const Eigen::Matrix<T, R, C>& matrix) {
      return matrix;
    }

    /**
     * Returns a matrix representation of the vector in column-major order
     * with the specified number of rows and columns.
     *
     * @tparam T type of the scalar
     * @param vec vector of values
     * @param m rows
     * @param n columns
     * @return the matrix representation of the input
     * @throw <code>std::invalid_argument</code> if the sizes do not match
     */
    template <typename T>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(const std::vector<T>& vec, int m, int n) {
      static const char* fun = "to_matrix(array)";
      check_size_match(fun, "rows * columns", m * n, "vector size", vec.size());
      return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic,
                                            Eigen::Dynamic> >(&vec[0], m, n);
    }

    /**
     * Returns a matrix representation of the standard vector of standard vectors
     * with the same dimensions and indexing order.
     *
     * @tparam T type of the scalar
     * @param vec vector of vectors of scalar values
     * @return the matrix representation of the input
     */
    template <typename T>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(const std::vector< std::vector<T> >& vec) {
      size_t R = vec.size();
      if (R != 0) {
        size_t C = vec[0].size();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(R, C);
        T* datap = result.data();
        for (size_t i=0, ij=0; i < C; i++)
          for (size_t j=0; j < R; j++, ij++)
            datap[ij] = vec[j][i];
        return result;
      } else {
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> (0, 0);
      }
    }


    // input is Eigen Row vector - Eigen::Matrix<T, 1, Eigen::Dynamic>
    // T is double, fvar, or var
    /**
     * Returns a matrix representation of the Eigen vector of standard vectors
     * with the same dimensions and indexing order.
     *
     * @tparam T type of the scalar
     * @param vec Eigen vector of vectors of scalar values
     * @return the matrix representation of the input
     */
    template <typename T>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(const std::vector< Eigen::Matrix<T, 1, Eigen::Dynamic> >& vec) {
      size_t R = vec.size();
      if (R != 0) {
        size_t C = vec[0].size();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(R, C);
        T* datap = result.data();
        for (size_t i=0, ij=0; i < C; i++)
          for (size_t j=0; j < R; j++, ij++)
            datap[ij] = vec[j][i];
        return result;
      } else {
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> (0, 0);
      }
    }

    /**
     * Returns a matrix representation of the standard vector of standard vectors
     * of integers with the same dimensions and indexing order.
     *
     * @param vec vector of vectors of integer values
     * @return the matrix representation of the input
     */
    inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(const std::vector< std::vector<int> >& vec) {
      size_t R = vec.size();
      if (R != 0) {
        size_t C = vec[0].size();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result(R, C);
        double* datap = result.data();
        for (size_t i=0, ij=0; i < C; i++)
          for (size_t j=0; j < R; j++, ij++)
            datap[ij] = vec[j][i];
        return result;
      } else {
        return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> (0, 0);
      }
    }

  }
}
#endif
