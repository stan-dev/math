#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_STD_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <test/unit/math/rev/mat/vectorize/build_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_autodiff.hpp>

template <typename F>
void expect_std_vector_value() {
  using std::vector;
  using stan::math::var;

  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    
    vector<var> y = build_vector<F>();
    vector<var> fy = F::template apply<vector<var> >(y);

    EXPECT_EQ(y.size(), fy.size());
    EXPECT_FLOAT_EQ(F::apply_base(y[i]).val(), fy[i].val());

    fy[i].grad();
    expect_autodiff<F>(y[i].val(), y[i].adj());
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {

    for (size_t j = 0; j < F::valid_inputs().size(); ++j) {

      vector<vector<var> > z;
      for (size_t i = 0; i < vector_vector_size; ++i) {
        z.push_back(build_vector<F>());
      } 
      vector<vector<var> > fz = 
                        F::template apply<vector<vector<var> > >(z);

      EXPECT_EQ(z.size(), fz.size());
      EXPECT_EQ(z[i].size(), fz[i].size());
      EXPECT_FLOAT_EQ(F::apply_base(z[i][j]).val(), fz[i][j].val());

      fz[i][j].grad();
      expect_autodiff<F>(z[i][j].val(), z[i][j].adj());
    }
  }
}    

#endif
