#include <stan/math/prim/prob/std_normal_lcdf.hpp>
#include <stan/math/prim/prob/std_normal_lpdf.hpp>
#include <gtest/gtest.h>
#include <array>

TEST(ProbStdNormal, scalar_tail_kernel) {
  // High-precision references: input, log Phi(input), and its slope.
  static constexpr std::array<double, 3> cases[]
      = {{-1e150, -4.9999999999999995e299, 1e150},
         {-30, -454.32124395634321, 30.033259667433676},
         {-6, -20.736768949974707, 6.1584826045445986},
         {-1, -1.8410216450092636, 1.5251352761609811},
         {0, -0.69314718055994529, 0.79788456080286541},
         {1, -0.17275377902344988, 0.28759997093917838},
         {5.6568542494923797, -7.7086289798515121e-09, 4.4895039573954219e-08},
         {5.6568542494923806, -7.7086289798514724e-09, 4.4895039573953994e-08},
         {6, -9.865876455243758e-10, 6.0758828558176762e-09},
         {30, -4.9067139271481872e-198, 1.4736461348785476e-196},
         {1e8, 0, 0},
         {1e150, 0, 0}};
  for (const auto& row : cases) {
    SCOPED_TRACE(row[0]);
    const auto result
        = stan::math::internal::std_normal_lcdf_value_grad<true>(row[0]);
    EXPECT_NEAR(row[1], result.first, 1e-12 * std::abs(row[1]));
    EXPECT_NEAR(row[2], result.second, 1e-12 * std::abs(row[2]));
    EXPECT_DOUBLE_EQ(
        result.first,
        stan::math::internal::std_normal_lcdf_value_grad<false>(row[0]).first);
  }
}

TEST(ProbStdNormal, large_finite_log_density) {
  EXPECT_TRUE(std::isfinite(stan::math::std_normal_lpdf(-1.5e154)));
  EXPECT_TRUE(std::isfinite(stan::math::std_normal_lcdf(-1.5e154)));
}

TEST(ProbStdNormal, vectorized_tail_kernel) {
  using stan::math::internal::std_normal_lcdf_value_grad;
  Eigen::ArrayXd z(23);
  z << -1e308, -1e150, -50, -6, -4 * stan::math::SQRT_TWO, -5, -1, -0.01, 0,
      0.01, 1, 3, 6, 20, 30, 37, 37.1, 38, 38.5, 40, 50, 1e150, 1e308;
  const auto result = std_normal_lcdf_value_grad<true>(z);
  const auto values = std_normal_lcdf_value_grad<false>(z);
  for (Eigen::Index i = 0; i < z.size(); ++i) {
    SCOPED_TRACE(z[i]);
    const auto scalar = std_normal_lcdf_value_grad<true>(z[i]);
    if (std::isfinite(scalar.first)) {
      EXPECT_NEAR(scalar.first, result.first[i],
                  std::abs(scalar.first) * 1e-12 + 1e-323);
    } else {
      EXPECT_EQ(scalar.first, result.first[i]);
    }
    EXPECT_EQ(result.first[i], values.first[i]);
    EXPECT_NEAR(scalar.second, result.second[i],
                std::abs(scalar.second) * 1e-12 + 1e-323);
  }
}

TEST(ProbStdNormal, integer_vectors) {
  using stan::math::std_normal_lcdf;
  const Eigen::Vector3i z(-1, 0, 1);
  EXPECT_NEAR(std_normal_lcdf(z.cast<double>().eval()), std_normal_lcdf(z),
              1e-14);
  EXPECT_NEAR(std_normal_lcdf(z), std_normal_lcdf(std::vector<int>{-1, 0, 1}),
              1e-14);
}
