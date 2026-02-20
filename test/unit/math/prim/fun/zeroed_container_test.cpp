#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/zeroed_container.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <type_traits>
#include <vector>

TEST(MathPrimFunZeroedContainer, autodiff_predicate_filters_and_zeros) {
  using stan::math::fvar;
  std::vector<fvar<double>> a{
      fvar<double>(1.3, 0.0),
      fvar<double>(-0.7, 0.0),
  };
  Eigen::Matrix<fvar<double>, -1, 1> b(2);
  b << fvar<double>(2.0, 0.0), fvar<double>(-3.0, 0.0);
  std::vector<int> x_i{5, 6};
  Eigen::VectorXd c{{4.0}};

  auto args = std::make_tuple(a, std::make_tuple(b, x_i), c);
  auto grads = stan::math::zeroed_container<stan::contains_autodiff>(args);

  using grads_t = std::decay_t<decltype(grads)>;
  EXPECT_EQ(2, std::tuple_size<grads_t>::value);

  const auto& a_grads = std::get<0>(grads);
  EXPECT_EQ(a.size(), a_grads.size());
  for (std::size_t i = 0; i < a_grads.size(); ++i) {
    EXPECT_FLOAT_EQ(0.0, a_grads[i]);
  }

  const auto& b_grads = std::get<0>(std::get<1>(grads));
  EXPECT_EQ(b.size(), b_grads.size());
  for (Eigen::Index i = 0; i < b_grads.size(); ++i) {
    EXPECT_FLOAT_EQ(0.0, b_grads(i));
  }
}

TEST(MathPrimFunZeroedContainer, integral_predicate_filters_and_zeros) {
  auto args = std::make_tuple(std::vector<int>{1, 2, 3}, 1.5,
                              std::make_tuple(2.5, 7));
  auto zeros = stan::math::zeroed_container<std::is_integral>(args);

  using zeros_t = std::decay_t<decltype(zeros)>;
  EXPECT_EQ(2, std::tuple_size<zeros_t>::value);

  const auto& vec_zero = std::get<0>(zeros);
  EXPECT_EQ(3, vec_zero.size());
  EXPECT_EQ(0, vec_zero[0]);
  EXPECT_EQ(0, vec_zero[1]);
  EXPECT_EQ(0, vec_zero[2]);

  const auto& nested_zero = std::get<0>(std::get<1>(zeros));
  EXPECT_EQ(0, nested_zero);
}
