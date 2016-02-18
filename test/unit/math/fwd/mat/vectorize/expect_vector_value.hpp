#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_VECTOR_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <vector>
#include <Eigen/Dense>
#include <test/unit/math/fwd/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fvar_eq.hpp>

template <typename F, typename T>
void expect_vector_value() {
  using stan::math::fvar;
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> vector_t;

  size_t num_inputs = F::valid_inputs().size();
  vector_t template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    vector_t b = build_matrix<F>(template_vector, i);
    vector_t fb = F::template apply<vector_t>(b);
    EXPECT_EQ(b.size(), fb.size());
    expect_fvar_eq(F::apply_base(b(i)), fb(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<vector_t> c;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i)
          c.push_back(build_matrix<F>(template_vector, j));
        else
          c.push_back(build_matrix<F>(template_vector));
      vector<vector_t> fc = F::template apply<vector<vector_t> >(c);

      EXPECT_EQ(c.size(), fc.size());
      EXPECT_EQ(c[i].size(), fc[i].size());
      expect_fvar_eq(F::apply_base(c[i](j)), fc[i](j));
    }
  }
}

#endif
