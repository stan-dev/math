#ifndef STAN_MATH_PRIM_MAT_FUN_TO_MATRIX_HPP
#define STAN_MATH_PRIM_MAT_FUN_TO_MATRIX_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
 // stan::scalar_type
#include <vector>

namespace stan {
  namespace math {

    // matrix to_matrix(matrix)
    // matrix to_matrix(vector)
    // matrix to_matrix(row_vector)
    template <typename T, int R, int C>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(Eigen::Matrix<T, R, C> matrix) {
      return matrix;
    }

    // matrix to_matrix(real[, ])
    template <typename T>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(const std::vector< std::vector<T> > & vec) {
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

    // matrix to_matrix(int[, ])
    inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    to_matrix(const std::vector< std::vector<int> > & vec) {
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
