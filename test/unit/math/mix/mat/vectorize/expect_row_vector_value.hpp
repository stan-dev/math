#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_ROW_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_ROW_VECTOR_VALUE_HPP

#include <vector>
#include <Eigen/Dense>
#include <test/unit/math/mix/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_fvar_var_eq.hpp>

template <typename F, typename T>
void expect_row_vector_value() {
  using std::vector;
  typedef Eigen::Matrix<T, 1, Eigen::Dynamic> row_vector_t;

  size_t num_inputs = F::valid_inputs().size();
  row_vector_t template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    row_vector_t a = build_matrix<F>(template_vector, i);
    row_vector_t b = build_matrix<F>(template_vector, i);
    row_vector_t fb = F::template apply<row_vector_t>(b);
    EXPECT_EQ(b.size(), fb.size());
    expect_fvar_var_eq(F::apply_base(a(i)), a(i), fb(i), b(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<row_vector_t> c;
      vector<row_vector_t> d;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i) { 
          c.push_back(build_matrix<F>(template_vector, j));
          d.push_back(build_matrix<F>(template_vector, j));
        }
        else {
          c.push_back(build_matrix<F>(template_vector));
          d.push_back(build_matrix<F>(template_vector));
        }
      vector<row_vector_t> fd = 
          F::template apply<vector<row_vector_t> >(d);

      EXPECT_EQ(d.size(), fd.size());
      EXPECT_EQ(d[i].size(), fd[i].size());
      expect_fvar_var_eq(
        F::apply_base(c[i](j)), c[i](j), fd[i](j), d[i](j));
    }
  }
}

#endif
