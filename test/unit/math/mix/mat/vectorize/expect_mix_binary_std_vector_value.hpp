#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_STD_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_std_vector_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
    bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
  for (size_t i = 0; i < template_v1.size(); ++i) {
    vector<input_t1> input_a1;
    vector<input_t1> input_a2;
    vector<input_t2> input_b1;
    vector<input_t2> input_b2;
    if (seed_one) {
      input_a1 = build_binary_vector1<F>(template_v1, i);
      input_a2 = build_binary_vector1<F>(template_v1, i);
    } else {
      input_a1 = build_binary_vector1<F>(template_v1);
      input_a2 = build_binary_vector1<F>(template_v1);
    }
    if (seed_two) {
      input_b1 = build_binary_vector2<F>(template_v2, i);
      input_b2 = build_binary_vector2<F>(template_v2, i);
    } else {
      input_b1 = build_binary_vector2<F>(template_v2);
      input_b2 = build_binary_vector2<F>(template_v2);
    }
    vector<FV> fa = F::template apply<vector<FV> >(input_a2, input_b2);
    EXPECT_EQ(input_a2.size(), fa.size());
    EXPECT_EQ(input_b2.size(), fa.size());
    expect_binary_val_deriv_eq(F::apply_base(input_a1[i], input_b1[i]),
                               input_a1[i], input_b1[i], fa[i], input_a2[i],
                               input_b2[i]);
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
    bool seed_one = 1, bool seed_two = 1) {
  using stan::math::var;
  using std::vector;

  const size_t num_v = 2;
  for (size_t i = 0; i < num_v; ++i) {
    for (size_t j = 0; j < template_v1.size(); ++j) {
      vector<vector<input_t1> > input_va1;
      vector<vector<input_t1> > input_va2;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i) {
          input_va1.push_back(build_binary_vector1<F>(template_v1, j));
          input_va2.push_back(build_binary_vector1<F>(template_v1, j));
        } else {
          input_va1.push_back(build_binary_vector1<F>(template_v1));
          input_va2.push_back(build_binary_vector1<F>(template_v1));
        }
      }
      vector<vector<input_t2> > input_vb1;
      vector<vector<input_t2> > input_vb2;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_two && k == i) {
          input_vb1.push_back(build_binary_vector2<F>(template_v2, j));
          input_vb2.push_back(build_binary_vector2<F>(template_v2, j));
        } else {
          input_vb1.push_back(build_binary_vector2<F>(template_v2));
          input_vb2.push_back(build_binary_vector2<F>(template_v2));
        }
      }
      vector<vector<FV> > fa
          = F::template apply<vector<vector<FV> > >(input_va2, input_vb2);
      EXPECT_EQ(input_va2.size(), fa.size());
      EXPECT_EQ(input_vb2.size(), fa.size());
      EXPECT_EQ(input_va2[i].size(), fa[i].size());
      EXPECT_EQ(input_vb2[i].size(), fa[i].size());
      expect_binary_val_deriv_eq(
          F::apply_base(input_va1[i][j], input_vb1[i][j]), input_va1[i][j],
          input_vb1[i][j], fa[i][j], input_va1[i][j], input_vb1[i][j]);
    }
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_std_vector_all_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(template_v1, template_v2);
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(template_v2, template_v1);
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_std_vector_std_vector_std_vector_all_eq(
    std::vector<input_t1> template_v1, std::vector<input_t2> template_v2) {
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      template_v1, template_v2);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      template_v2, template_v1);
}

template <typename F, typename FV>
void expect_mix_binary_std_vector_value() {
  using stan::math::fvar;
  using std::vector;

  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<FV> var_template_v;

  expect_mix_binary_std_vector_std_vector_all_eq<F, FV>(var_template_v,
                                                        int_template_v);
  expect_mix_binary_std_vector_std_vector_all_eq<F, FV>(var_template_v,
                                                        d_template_v);
  expect_mix_binary_std_vector_std_vector_all_eq<F, FV>(var_template_v,
                                                        var_template_v);
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(var_template_v,
                                                    var_template_v, 1, 0);
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(var_template_v,
                                                    var_template_v, 0, 1);

  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_all_eq<F, FV>(
      var_template_v, int_template_v);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_all_eq<F, FV>(
      var_template_v, d_template_v);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_all_eq<F, FV>(
      var_template_v, var_template_v);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      var_template_v, var_template_v, 1, 0);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
      var_template_v, var_template_v, 0, 1);
}
#endif
