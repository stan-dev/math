#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BUILD_TEMPLATE_MATRIX_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BUILD_TEMPLATE_MATRIX_HPP

#include <Eigen/Dense>

template <typename matrix_t> static inline
Eigen::Matrix<matrix_t, Eigen::Dynamic, Eigen::Dynamic>
build_template_matrix(Eigen::Matrix<matrix_t, Eigen::Dynamic,
                                    Eigen::Dynamic> template_m,
                      int size_1, int size_2) {
  template_m.resize(size_1, size_2);
  return template_m;
}

template <typename matrix_t> static inline
Eigen::Matrix<matrix_t, 1, Eigen::Dynamic> build_template_matrix(
    Eigen::Matrix<matrix_t, 1, Eigen::Dynamic> template_m,
    int size_1, int size_2) {
  template_m.resize(size_1);
  return template_m;
}

template <typename matrix_t> static inline
Eigen::Matrix<matrix_t, Eigen::Dynamic, 1> build_template_matrix(
    Eigen::Matrix<matrix_t, Eigen::Dynamic, 1> template_m,
    int size_1, int size_2) {
  template_m.resize(size_1);
  return template_m;
}

#endif
