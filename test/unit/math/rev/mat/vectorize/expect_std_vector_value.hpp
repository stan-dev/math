#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_STD_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <test/unit/math/rev/mat/vectorize/build_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_eq.hpp>

template <typename F>
void expect_std_vector_value() {
  using std::vector;
  using stan::math::var;

  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    
    vector<var> y = build_vector<F>();
    vector<var> fy = F::template apply<vector<var> >(y);
    vector<var> z = build_vector<F>();

    EXPECT_EQ(y.size(), fy.size());
    expect_eq(F::apply_base(z[i]), z[i], fy[i], y[i]);
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {

    for (size_t j = 0; j < F::valid_inputs().size(); ++j) {

      vector<vector<var> > a;
      for (size_t i = 0; i < vector_vector_size; ++i) {
        a.push_back(build_vector<F>());
      } 
      vector<vector<var> > fa = 
                        F::template apply<vector<vector<var> > >(a);
      EXPECT_EQ(a.size(), fa.size());
      EXPECT_EQ(a[i].size(), fa[i].size());

      vector<vector<var> > b;
      for (size_t i = 0; i < vector_vector_size; ++i)
        b.push_back(build_vector<F>());

      expect_eq(F::apply_base(b[i][j]), b[i][j], fa[i][j], a[i][j]);
    }
  }
}    

#endif
