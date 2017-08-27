#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_STD_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_binary_val_deriv_eq.hpp>

template <typename F, typename input_t1, typename input_t2> void
expect_rev_binary_std_vector_std_vector_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {
  using std::vector;
  using stan::math::var;

  for (size_t i = 0; i < template_v1.size(); ++i) {
    vector<input_t1> input_va1 = build_binary_vector1<F>(template_v1);
    vector<input_t1> input_va2 = build_binary_vector1<F>(template_v1);
    vector<input_t2> input_vb1 = build_binary_vector2<F>(template_v2);
    vector<input_t2> input_vb2 = build_binary_vector2<F>(template_v2);
    vector<var> fa = F::template apply<vector<var> >(input_va2, input_vb2);
    EXPECT_EQ(input_va2.size(), fa.size());
    EXPECT_EQ(input_vb2.size(), fa.size());
    expect_binary_val_deriv_eq(F::apply_base(input_va1[i], input_vb1[i]),
                               input_va1[i], input_vb1[i],
                               fa[i], input_va2[i], input_vb2[i]);
  }
}

template <typename F, typename input_t1, typename input_t2> void
expect_rev_binary_std_vector_std_vector_std_vector_std_vector_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {

  using std::vector;
  using stan::math::var;

  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < template_v1.size(); ++j) {
      vector<vector<input_t1> > input_va1;
      input_va1.push_back(build_binary_vector1<F>(template_v1));
      input_va1.push_back(build_binary_vector1<F>(template_v1));
      vector<vector<input_t1> > input_va2;
      input_va2.push_back(build_binary_vector1<F>(template_v1));
      input_va2.push_back(build_binary_vector1<F>(template_v1));
      vector<vector<input_t2> > input_vb1;
      input_vb1.push_back(build_binary_vector2<F>(template_v2));
      input_vb1.push_back(build_binary_vector2<F>(template_v2));
      vector<vector<input_t2> > input_vb2;
      input_vb2.push_back(build_binary_vector2<F>(template_v2));
      input_vb2.push_back(build_binary_vector2<F>(template_v2));
      vector<vector<var> > fa = F::template apply<vector<vector<var> > >(
        input_va2, input_vb2);
      EXPECT_EQ(input_va2.size(), fa.size());
      EXPECT_EQ(input_vb2.size(), fa.size());
      EXPECT_EQ(input_va2[i].size(), fa[i].size());
      EXPECT_EQ(input_vb2[i].size(), fa[i].size());
      expect_binary_val_deriv_eq(F::apply_base(input_va1[i][j],
                                               input_vb1[i][j]),
                                 input_va1[i][j], input_vb1[i][j], fa[i][j],
                                 input_va2[i][j], input_vb2[i][j]);
    }
  }
}

template <typename F, typename input_t1, typename input_t2> void
expect_rev_binary_std_vector_std_vector_all_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {
  expect_rev_binary_std_vector_std_vector_eq<F>(template_v1, template_v2);
  expect_rev_binary_std_vector_std_vector_eq<F>(template_v2, template_v1);
}

template <typename F, typename input_t1, typename input_t2> void
expect_rev_binary_std_vector_std_vector_std_vector_std_vector_all_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {
  expect_rev_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(
      template_v1, template_v2);
  expect_rev_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(
      template_v2, template_v1);
}

template <typename F>
void expect_rev_binary_std_vector_value() {
  using std::vector;
  using stan::math::var;

  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<var> var_template_v;

  expect_rev_binary_std_vector_std_vector_all_eq<F>(var_template_v,
                                                    int_template_v);
  expect_rev_binary_std_vector_std_vector_all_eq<F>(var_template_v,
                                                    d_template_v);
  expect_rev_binary_std_vector_std_vector_all_eq<F>(var_template_v,
                                                    var_template_v);

  expect_rev_binary_std_vector_std_vector_std_vector_std_vector_all_eq<F>(
      var_template_v, int_template_v);
  expect_rev_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(
      var_template_v, d_template_v);
  expect_rev_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(
      var_template_v, var_template_v);
}
#endif
