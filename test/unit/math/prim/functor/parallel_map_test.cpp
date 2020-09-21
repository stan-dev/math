#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <gtest/gtest.h>

TEST(MathFunctions, parall_map_ranged) {
  using Eigen::VectorXd;
  Eigen::VectorXd in1_par = Eigen::VectorXd::Random(1000);
  Eigen::VectorXd in1_ser = in1_par;
  Eigen::VectorXd out_par(1000);
  Eigen::VectorXd out_ser(1000);


  // Functor defining how inputs should be indexed
  auto ind_f = [&](int begin, int end, const auto& fun,
                    const auto& x) {
    return fun(x.segment(begin, end-begin));
  };

  // Functor defining function to be applied to indexed arguments
  auto f = [&](const auto& x) {
    return x.array().exp().matrix(); };

  // Boolean template parameter to enable ranged parallelism
  stan::math::parallel_map<true>(f, ind_f, std::forward<Eigen::VectorXd>(out_par),
                         std::forward_as_tuple(in1_par));
  EXPECT_MATRIX_FLOAT_EQ(out_par, in1_par.array().exp().matrix());
}

