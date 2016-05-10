#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_SCALAR_EQ
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_SCALAR_EQ

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_binary_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename FV, typename input_t1, typename input_t2, 
int R, int C> void expect_fwd_binary_matrix_scalar_eq(const Eigen::Matrix<
input_t1, R, C>& template_m, const std::vector<input_t2>& template_scalar_v,
bool seed_one = 1, bool seed_two = 1) {
  using std::vector;

  typedef Eigen::Matrix<input_t1, R, C> input_matrix_t;
  typedef Eigen::Matrix<FV, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    for (int j = 0; j < template_m.size(); ++j) {
      input_matrix_t input_m; 
      if (seed_one)
        input_m = build_fwd_binary_matrix1<F>(i, template_m, j);
      else
        input_m = build_fwd_binary_matrix1<F>(i, template_m);
      input_t2 input_2;
      if (seed_two)
        input_2 = build_binary_vector2<F>(template_scalar_v, i)[i];
      else
        input_2 = build_binary_vector2<F>(template_scalar_v)[i];
      result_matrix_t fa = F::template apply<result_matrix_t>(
      input_m, input_2);
      EXPECT_EQ(input_m.size(), fa.size());
      expect_val_deriv_eq(F::apply_base(input_m(j), input_2), fa(j));
    } 
  }   
}
#endif
