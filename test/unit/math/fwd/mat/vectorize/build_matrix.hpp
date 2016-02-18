#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_MATRIX_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_MATRIX_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <Eigen/Dense>
#include <test/unit/math/fwd/mat/vectorize/build_vector.hpp>

template <typename F, typename T>
static inline Eigen::Matrix<typename Eigen::internal::traits<T>::Scalar, 
                              T::RowsAtCompileTime,
                                T::ColsAtCompileTime>
build_matrix(const T& x, int seed_index = -1) {

  typedef typename Eigen::internal::traits<T>::Scalar fvar_type;
  Eigen::Matrix<fvar_type, T::RowsAtCompileTime, T::ColsAtCompileTime>
    fvar_matrix(x.rows(), x.cols());
  size_t num_inputs = F::valid_inputs().size();
  for (int i = 0; i < x.size(); ++i) {
    std::vector<fvar_type> inputs;
    if (seed_index == i)
      inputs = build_vector<F>(std::vector<fvar_type>(), 
                                       (seed_index % num_inputs));
    else
      inputs = build_vector<F>(std::vector<fvar_type>()); 
    fvar_matrix(i) = inputs[(i % num_inputs)];
  }
  return fvar_matrix;
}

#endif
