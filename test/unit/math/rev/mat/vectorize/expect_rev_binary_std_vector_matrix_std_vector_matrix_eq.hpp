#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/build_rev_binary_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename matrix_t1, typename matrix_t2>
void expect_rev_binary_std_vector_matrix_std_vector_matrix_eq(
const matrix_t1& template_m1, const matrix_t2& template_m2) {
  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<var, matrix_t1::RowsAtCompileTime, 
  matrix_t1::ColsAtCompileTime> result_matrix_t;

  for (size_t i = 0; i < 2; ++i) {
    for (int j = 0; j < template_m1.size(); ++j) {
      vector<matrix_t1> input_mva1;
      input_mva1.push_back(build_rev_binary_matrix1<F>(template_m1));
      input_mva1.push_back(build_rev_binary_matrix1<F>(template_m1));
      vector<matrix_t1> input_mva2;
      input_mva2.push_back(build_rev_binary_matrix1<F>(template_m1));
      input_mva2.push_back(build_rev_binary_matrix1<F>(template_m1));
      vector<matrix_t2> input_mvb1;
      input_mvb1.push_back(build_rev_binary_matrix2<F>(template_m2));
      input_mvb1.push_back(build_rev_binary_matrix2<F>(template_m2));
      vector<matrix_t2> input_mvb2;
      input_mvb2.push_back(build_rev_binary_matrix2<F>(template_m2));
      input_mvb2.push_back(build_rev_binary_matrix2<F>(template_m2));
      vector<result_matrix_t> fa = F::template 
      apply<vector<result_matrix_t> >(input_mva2, input_mvb2);
      EXPECT_EQ(input_mva2.size(), fa.size());
      EXPECT_EQ(input_mvb2.size(), fa.size());
      EXPECT_EQ(input_mva2[i].size(), fa[i].size());
      EXPECT_EQ(input_mvb2[i].size(), fa[i].size());
      expect_binary_val_deriv_eq(F::apply_base(input_mva1[i](j), 
      input_mvb1[i](j)), input_mva1[i](j), input_mvb1[i](j), 
      fa[i](j), input_mva2[i](j), input_mvb2[i](j));
    }
  }   
}
#endif
