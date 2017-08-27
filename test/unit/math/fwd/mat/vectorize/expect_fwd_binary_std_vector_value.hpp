#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_STD_VECTOR_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <vector>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_fwd_binary_std_vector_std_vector_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
    bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
  for (size_t i = 0; i < template_v1.size(); ++i) {
    vector<input_t1> input_a;
    vector<input_t2> input_b;
    if (seed_one)
      input_a = build_binary_vector1<F>(template_v1, i);
    else
      input_a = build_binary_vector1<F>(template_v1);
    if (seed_two)
      input_b = build_binary_vector2<F>(template_v2, i);
    else
      input_b = build_binary_vector2<F>(template_v2);
    vector<FV> fa = F::template apply<vector<FV> >(input_a, input_b);
    EXPECT_EQ(input_a.size(), fa.size());
    EXPECT_EQ(input_b.size(), fa.size());
    expect_val_deriv_eq(F::apply_base(input_a[i], input_b[i]), fa[i]);
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
    bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
  using stan::math::var;

  const size_t num_v = 2;
  for (size_t i = 0; i < num_v; ++i) {
    for (size_t j = 0; j < template_v1.size(); ++j) {
      vector<vector<input_t1> > input_v1;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i)
          input_v1.push_back(build_binary_vector1<F>(template_v1, j));
        else
          input_v1.push_back(build_binary_vector1<F>(template_v1));
      }
      vector<vector<input_t2> > input_v2;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_two && k == i)
          input_v2.push_back(build_binary_vector2<F>(template_v2, j));
        else
          input_v2.push_back(build_binary_vector2<F>(template_v2));
      }
      vector<vector<FV> > fa = F::template apply<vector<vector<FV> > >(
          input_v1, input_v2);
      EXPECT_EQ(input_v1.size(), fa.size());
      EXPECT_EQ(input_v2.size(), fa.size());
      EXPECT_EQ(input_v1[i].size(), fa[i].size());
      EXPECT_EQ(input_v2[i].size(), fa[i].size());
      expect_val_deriv_eq(F::apply_base(input_v1[i][j], input_v2[i][j]),
                          fa[i][j]);
    }
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_fwd_binary_std_vector_std_vector_all_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {
  expect_fwd_binary_std_vector_std_vector_eq<F, FV>(template_v1,
                                                    template_v2);
  expect_fwd_binary_std_vector_std_vector_eq<F, FV>(template_v2,
                                                    template_v1);
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_all_eq(
std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {
  expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      template_v1, template_v2);
  expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      template_v2, template_v1);
}

template <typename F, typename FV>
void expect_fwd_binary_std_vector_value() {
  using std::vector;
  using stan::math::fvar;

  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<FV> var_template_v;

  expect_fwd_binary_std_vector_std_vector_all_eq<F, FV>(var_template_v,
                                                        int_template_v);
  expect_fwd_binary_std_vector_std_vector_all_eq<F, FV>(var_template_v,
                                                        d_template_v);
  expect_fwd_binary_std_vector_std_vector_all_eq<F, FV>(var_template_v,
                                                        var_template_v);
  expect_fwd_binary_std_vector_std_vector_eq<F, FV>(var_template_v,
                                                    var_template_v, 1, 0);
  expect_fwd_binary_std_vector_std_vector_eq<F, FV>(var_template_v,
                                                    var_template_v, 0, 1);

  expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_all_eq
      <F, FV>(var_template_v, int_template_v);
  expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_all_eq
      <F, FV>(var_template_v, d_template_v);
  expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_all_eq
      <F, FV>(var_template_v, var_template_v);
  expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      var_template_v, var_template_v, 1, 0);
  expect_fwd_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      var_template_v, var_template_v, 0, 1);
}
#endif
