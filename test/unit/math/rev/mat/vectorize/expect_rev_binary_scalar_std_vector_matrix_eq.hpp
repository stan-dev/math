#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_SCALAR_STD_VECTOR_MATRIX_EQ
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_SCALAR_STD_VECTOR_MATRIX_EQ

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename input_t1, typename input_t2, int R, int C>
void expect_rev_binary_scalar_std_vector_matrix_eq(
const std::vector<input_t1>& template_scalar_v,
const Eigen::Matrix<input_t2, R, C>& template_m) {
  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<input_t2, R, C> input_matrix_t;
  typedef Eigen::Matrix<var, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < template_m.size(); ++k) {
        vector<input_t1> input_va1 = build_binary_vector1<F>(
        template_scalar_v);
        vector<input_t1> input_va2 = build_binary_vector1<F>(
        template_scalar_v);
        vector<input_t2> input_vb1 = build_binary_vector2<F>(
        vector<input_t2>());
        vector<input_t2> input_vb2 = build_binary_vector2<F>(
        vector<input_t2>());
        vector<input_matrix_t> input_mv1;
        input_mv1.push_back(build_binary_matrix(input_vb1[i], template_m));
        input_mv1.push_back(build_binary_matrix(input_vb1[i], template_m));
        vector<input_matrix_t> input_mv2;
        input_mv2.push_back(build_binary_matrix(input_vb2[i], template_m));
        input_mv2.push_back(build_binary_matrix(input_vb2[i], template_m));
        vector<result_matrix_t> fa = F::template 
        apply<vector<result_matrix_t> >(input_va2[i], input_mv2);
        EXPECT_EQ(input_mv2.size(), fa.size());
        EXPECT_EQ(input_mv2[j].size(), fa[j].size());
        expect_binary_val_deriv_eq(F::apply_base(input_va1[i], 
        input_mv1[j](k)), input_va1[i], input_mv1[j](k), 
        fa[j](k), input_va2[i], input_mv2[j](k));
      } 
    }   
  }
}
#endif
