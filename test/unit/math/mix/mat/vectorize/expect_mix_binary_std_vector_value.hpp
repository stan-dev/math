#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_BINARY_STD_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_scalar_std_vector_eq(
std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
//  for (size_t i = 0; i < template_v1.size(); ++i) {
//    for (size_t j = 0; j < 5; ++j) {
      size_t i = 0;
      size_t j = 0;
      input_t1 input_a1 = 0;
      input_t1 input_a2 = 0;
      if (seed_one) { 
        input_a1 = build_binary_vector1<F>(template_v1, i)[i];
        input_a2 = build_binary_vector1<F>(template_v1, i)[i];
      } 
      else {
        input_a1 = build_binary_vector1<F>(template_v1)[i];
        input_a2 = build_binary_vector1<F>(template_v1)[i];
      }
      vector<input_t2> input_v1;
      vector<input_t2> input_v2;
      for (size_t k = 0; k < 5; ++k) {
        input_v1.push_back(build_binary_vector2<F>(template_v2)[i]);
        input_v2.push_back(build_binary_vector2<F>(template_v2)[i]);
      }
      if (seed_two) {
        input_v1[j] = build_binary_vector2<F>(template_v2, i)[i]; 
        input_v2[j] = build_binary_vector2<F>(template_v2, i)[i]; 
      }
      vector<FV> fa = F::template apply<vector<FV> >(input_a2, input_v2);
      EXPECT_EQ(input_v2.size(), fa.size());
      //std::cout << "Well" << j << std::endl;
      //std::cout << input_a1 << std::endl;
      //std::cout << input_v1[j] << std::endl;
      expect_binary_val_deriv_eq(F::apply_base(input_a1, input_v1[j]), 
      input_a1, input_v1[j], fa[j], input_a2, input_v2[j]);
      expect_binary_val_deriv_eq(F::apply_base(input_a1, input_v1[j]), 
      input_a1, input_v1[j], fa[j], input_a2, input_v2[j]);
     /*
      FV fb = F::apply_base(input_a1, input_v1[j]);
      AVEC exp_b = createAVEC(input_a1.d_);
      VEC exp_bg;
      fb.d_.grad(exp_b, exp_bg);
      std::cout << input_a1 << std::endl;
      std::cout << input_v1[j] << std::endl;
      std::cout << exp_bg[0] << std::endl;
      stan::math::set_zero_all_adjoints();
      std::cout << exp_bg[0] << std::endl;
      fa[j].d_.grad();
      std::cout << exp_bg[0] << std::endl;
  std::cout << input_a1.val() << std::endl;
      std::cout << input_a2.val() << std::endl;
      std::cout << input_v1[j] << std::endl;
      std::cout << input_v2[j] << std::endl;
      //F::template apply<vector<FV> >(input_a1, input_v1);
      std::cout << fa[j] << std::endl;
      std::cout << fa[j].d_ << std::endl;
      FV ab = F::apply_base(input_a1, input_v1[j]);
      std::cout << ab.val() << std::endl;
      std::cout << ab.d_ << std::endl;
      expect_binary_val_deriv_eq(F::apply_base(input_a1, input_v1[j]), 
      input_a1, input_v1[j], fa[j], input_a2, input_v2[j]);
      std::cout << input_a1.val() << std::endl;
      std::cout << input_v1[j] << std::endl;
      //FV ab = F::apply_base(input_a1, input_v1[j]);
      expect_binary_val_deriv_eq(F::apply_base(input_a1, input_v1[j]), 
      input_a1, input_v1[j], fa[j], input_a2, input_v2[j]);
      std::cout << "Yay" << std::endl;
      expect_binary_val_deriv_eq(F::apply_base(input_a1, input_v1[j]), 
      input_a1, input_v1[j], fa[j], input_a2, input_v2[j]);
*/
//    }
//  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_scalar_eq(
std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
  for (size_t i = 0; i < template_v2.size(); ++i) {
    for (size_t j = 0; j < 5; ++j) {
      vector<input_t1> input_v1(5, build_binary_vector1<F>(template_v1)[i]);
      vector<input_t1> input_v2(5, build_binary_vector1<F>(template_v1)[i]);
      if (seed_one) {
        input_v1[j] = build_binary_vector1<F>(template_v1, i)[i];
        input_v2[j] = build_binary_vector1<F>(template_v1, i)[i];
      }
      input_t2 input_b1 = 0; 
      input_t2 input_b2 = 0; 
      if (seed_two) {
        input_b1 = build_binary_vector2<F>(template_v2, i)[i];
        input_b2 = build_binary_vector2<F>(template_v2, i)[i];
      }
      else {
        input_b1 = build_binary_vector2<F>(template_v2)[i];
        input_b2 = build_binary_vector2<F>(template_v2)[i];
      }
      vector<FV> fa = F::template apply<vector<FV> >(input_v2, input_b2);
      EXPECT_EQ(input_v2.size(), fa.size());
      expect_binary_val_deriv_eq(F::apply_base(input_v1[j], input_b1),
      input_v1[j], input_b1, fa[j], input_v2[j], input_b2);
    }
  }
}
/*
template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_scalar_std_vector_std_vector_eq(
std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
  using stan::math::var;

  const size_t num_v = 2;
  for (size_t i = 0; i < template_v1.size(); ++i) {
    for (size_t j = 0; j < num_v; ++j) {
      for (size_t k = 0; k < 5; ++k) {
        input_t1 input_a;
        if (seed_one)
          input_a = build_binary_vector1<F>(template_v1, i)[i];
        else
          input_a = build_binary_vector1<F>(template_v1)[i];
        vector<vector<input_t2> > input_v;
        for (size_t l = 0; l < num_v; ++l) {
          vector<input_t2> input_b(5, build_binary_vector2<F>(
          template_v2)[i]);
          if (seed_two && l == j)
            input_b[k] = build_binary_vector2<F>(template_v2, i)[i];
          input_v.push_back(input_b);
        }
        vector<vector<FV> > fa = F::template apply<vector<vector<FV> > >(
        input_a, input_v);
        EXPECT_EQ(input_v.size(), fa.size());
        EXPECT_EQ(input_v[j].size(), fa[j].size());
        expect_binary_val_deriv_eq(F::apply_base(input_a, input_v[j][k]), 
        fa[j][k]);
      }
    }
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_std_vector_scalar_eq(
std::vector<input_t1> template_v1, std::vector<input_t2> template_v2,
bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
  using stan::math::var;

  const size_t num_v = 2;
  for (size_t i = 0; i < template_v2.size(); ++i) {
    for (size_t j = 0; j < num_v; ++j) {
      for (size_t k = 0; k < 5; ++k) {
        vector<vector<input_t1> > input_v;
        for (size_t l = 0; l < num_v; ++l) { 
          vector<input_t1> input_a(5, build_binary_vector2<F>(
          template_v1)[i]);
          if (seed_one && l == j)
            input_a[k] = build_binary_vector2<F>(template_v1, i)[i];
          input_v.push_back(input_a);
        }
        input_t2 input_b;
        if (seed_two)
          input_b = build_binary_vector2<F>(template_v2, i)[i];
        else
          input_b = build_binary_vector2<F>(template_v2)[i];
        vector<vector<FV> > fa = F::template apply<vector<vector<FV> > >(
        input_v, input_b);
        EXPECT_EQ(input_v.size(), fa.size());
        EXPECT_EQ(input_v[j].size(), fa[j].size());
        expect_binary_val_deriv_eq(F::apply_base(input_v[j][k], input_b), 
        fa[j][k]);
      }
    }
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_std_vector_eq(
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
    expect_binary_val_deriv_eq(F::apply_base(input_a[i], input_b[i]), fa[i]);
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2>
void expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq(
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
      expect_binary_val_deriv_eq(F::apply_base(input_v1[i][j], input_v2[i][j]), 
      fa[i][j]);
    }
  }
}
*/
template <typename F, typename FV>
void expect_mix_binary_std_vector_value() {
  using std::vector;
  using stan::math::fvar;

  vector<int> int_template_v(F::int_valid_inputs1().size());
  vector<double> d_template_v(F::valid_inputs1().size());
  vector<FV> var_template_v(F::valid_inputs1().size());
/*
  for (size_t i = 0; i < template_v1.size(); ++i) {
    for (size_t j = 0; j < 5; ++j) {
      input_t1 input_a1 = 0;
      input_t1 input_a2 = 0;
      if (seed_one) { 
        input_a1 = build_binary_vector1<F>(template_v1, i)[i];
        input_a2 = build_binary_vector1<F>(template_v1, i)[i];
      } 
      else {
        input_a1 = build_binary_vector1<F>(template_v1)[i];
        input_a2 = build_binary_vector1<F>(template_v1)[i];
      }
      vector<input_t2> input_v1(5, build_binary_vector2<F>(template_v2)[i]);
      vector<input_t2> input_v2(5, build_binary_vector2<F>(template_v2)[i]);
      if (seed_two) {
        input_v1[j] = build_binary_vector2<F>(template_v2, i)[i]; 
        input_v2[j] = build_binary_vector2<F>(template_v2, i)[i]; 
      }
      vector<FV> fa = F::template apply<vector<FV> >(input_a2, input_v2);
      EXPECT_EQ(input_v2.size(), fa.size());
      expect_binary_val_deriv_eq(F::apply_base(input_a1, input_v1[j]), 
      input_a1, input_v1[j], fa[j], input_a2, input_v2[j]);
    }
  } 
  //scalar, vector
  //var, int
*/
  expect_mix_binary_scalar_std_vector_eq<F, FV>(var_template_v, 
  int_template_v);
/*
  expect_mix_binary_std_vector_scalar_eq<F, FV>(var_template_v, 
  int_template_v);
  expect_mix_binary_scalar_std_vector_eq<F, FV>(int_template_v, 
  var_template_v);
  expect_mix_binary_std_vector_scalar_eq<F, FV>(int_template_v, 
  var_template_v);
  //var, double
  expect_mix_binary_scalar_std_vector_eq<F, FV>(var_template_v, 
  d_template_v);
  expect_mix_binary_std_vector_scalar_eq<F, FV>(var_template_v, 
  d_template_v);
  expect_mix_binary_scalar_std_vector_eq<F, FV>(d_template_v, 
  var_template_v);
  expect_mix_binary_std_vector_scalar_eq<F, FV>(d_template_v, 
  var_template_v);
  //var, var
  expect_mix_binary_scalar_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 1, 0);
  expect_mix_binary_scalar_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 0, 1);
  expect_mix_binary_scalar_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 1, 1);
  expect_mix_binary_std_vector_scalar_eq<F, FV>(var_template_v, 
  var_template_v, 1, 0);
  expect_mix_binary_std_vector_scalar_eq<F, FV>(var_template_v, 
  var_template_v, 0, 1);
  expect_mix_binary_std_vector_scalar_eq<F, FV>(var_template_v, 
  var_template_v, 1, 1);
  
  //vector, vector
  //var, int
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(var_template_v, 
  int_template_v);
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(int_template_v, 
  var_template_v);
  //var, double
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(var_template_v, 
  d_template_v);
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(d_template_v, 
  var_template_v);
  //var, var
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 1, 0);
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 0, 1);
  expect_mix_binary_std_vector_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 1, 1);

  //scalar, vector<vector>
  //var, int 
  expect_mix_binary_scalar_std_vector_std_vector_eq<F, FV>(var_template_v, 
  int_template_v);
  expect_mix_binary_std_vector_std_vector_scalar_eq<F, FV>(var_template_v, 
  int_template_v);
  expect_mix_binary_scalar_std_vector_std_vector_eq<F, FV>(int_template_v, 
  var_template_v);
  expect_mix_binary_std_vector_std_vector_scalar_eq<F, FV>(int_template_v, 
  var_template_v);
  //var, double 
  expect_mix_binary_scalar_std_vector_std_vector_eq<F, FV>(var_template_v, 
  d_template_v);
  expect_mix_binary_std_vector_std_vector_scalar_eq<F, FV>(var_template_v, 
  d_template_v);
  expect_mix_binary_scalar_std_vector_std_vector_eq<F, FV>(d_template_v, 
  var_template_v);
  expect_mix_binary_std_vector_std_vector_scalar_eq<F, FV>(d_template_v, 
  var_template_v);
  //var, var 
  expect_mix_binary_scalar_std_vector_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 1, 0);
  expect_mix_binary_scalar_std_vector_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 0, 1);
  expect_mix_binary_scalar_std_vector_std_vector_eq<F, FV>(var_template_v, 
  var_template_v, 1, 1);
  expect_mix_binary_std_vector_std_vector_scalar_eq<F, FV>(var_template_v, 
  var_template_v, 0, 1);
  expect_mix_binary_std_vector_std_vector_scalar_eq<F, FV>(var_template_v, 
  var_template_v, 1, 0);
  expect_mix_binary_std_vector_std_vector_scalar_eq<F, FV>(var_template_v, 
  var_template_v, 1, 1);
 
  //vector<vector>, vector<vector>
  //var, int
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
  var_template_v, int_template_v);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
  int_template_v, var_template_v);
  //var, double
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
  var_template_v, d_template_v);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
  d_template_v, var_template_v);
  //var, var 
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
  var_template_v, var_template_v, 1, 0);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
  var_template_v, var_template_v, 0, 1);
  expect_mix_binary_std_vector_std_vector_std_vector_std_vector_eq<F, FV>(
  var_template_v, var_template_v, 1, 1);
*/
}    
#endif
