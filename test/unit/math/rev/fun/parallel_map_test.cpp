#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <gtest/gtest.h>

TEST(MathFunctions, parall_map) {
  using stan::math::pow;
  using stan::math::var;
  using stan::math::vector_v;
  vector_v in1_par = vector_v::Random(1000);
  vector_v in2_par = vector_v::Random(1000);
  vector_v in1_ser = in1_par;
  vector_v in2_ser = in2_par;
  vector_v out_par(1000);
  vector_v out_ser(1000);
  
  for(int i = 0; i < 1000; ++i) {
    out_ser[i] = in1_ser[i] * 0.5 + exp(in2_ser[i]);
  }

  // Functor defining how inputs should be indexed
  auto ind_f = [&](int i, const auto& fun,
                    const auto& x, const auto& y, const auto& z) {
    return fun(x.coeffRef(i), y, z.coeffRef(i));
  };

  // Functor defining function to be applied to indexed arguments
  auto f = [&](const auto& x, const auto& y, const auto& z) {
    return x * y + exp(z); };

  out_par = parallel_map(f, ind_f, std::forward<vector_v>(out_par),
                         std::forward_as_tuple(in1_par,0.5,in2_par));
  EXPECT_MATRIX_EQ(out_par.val(), out_ser.val());
  EXPECT_MATRIX_EQ(out_par.adj(), out_ser.adj());

  for(int i = 0; i < 1000; ++i) {
    out_par[i].grad();
    out_ser[i].grad();
  }

  EXPECT_MATRIX_EQ(in1_par.adj(), in1_ser.adj());
  EXPECT_MATRIX_EQ(in2_par.adj(), in2_ser.adj());
}

TEST(MathFunctions, parall_map_vec) {
  using stan::math::pow;
  using stan::math::var;
  using stan::math::vector_v;
  vector_v in1_par = vector_v::Random(1000);
  vector_v in2_par = vector_v::Random(1000);
  vector_v in3_par = vector_v::Random(1000);
  vector_v in1_ser = in1_par;
  vector_v in2_ser = in2_par;
  vector_v in3_ser = in3_par;
  vector_v out_par(1000);
  vector_v out_ser(1000);
  
  for(int i = 0; i < 1000; ++i) {
    out_ser[i] = in1_ser[i] * sum(in2_ser) + exp(in3_ser[i]);
  }

  // Functor defining how inputs should be indexed
  auto ind_f = [&](int i, const auto& fun,
                    const auto& x, const auto& y, const auto& z) {
    return fun(x.coeffRef(i), y, z.coeffRef(i));
  };

  // Functor defining function to be applied to indexed arguments
  auto f = [&](const auto& x, const auto& y, const auto& z) {
    return x * sum(y) + exp(z); };

  out_par = parallel_map(f, ind_f, std::forward<vector_v>(out_par),
                         std::forward_as_tuple(in1_par,in2_par,in3_par));
  EXPECT_MATRIX_EQ(out_par.val(), out_ser.val());
  EXPECT_MATRIX_EQ(out_par.adj(), out_ser.adj());

  for(int i = 0; i < 1000; ++i) {
    out_par[i].grad();
    out_ser[i].grad();
  }

  EXPECT_MATRIX_EQ(in1_par.adj(), in1_ser.adj());
  EXPECT_MATRIX_EQ(in2_par.adj(), in2_ser.adj());
  EXPECT_MATRIX_EQ(in3_par.adj(), in3_ser.adj());
}
