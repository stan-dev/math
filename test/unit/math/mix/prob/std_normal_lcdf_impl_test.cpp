#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/prob/normal_lcdf_tail_test_helpers.hpp>
#include <array>
#include <limits>

using normal_lcdf_tail_test::check_tail_derivatives;

TEST_F(AgradRev, std_normal_extreme_tail_derivatives) {
  using namespace stan::math;
  // High-precision references for L', L'', L''', where L(z)=log Phi(z).
  // The last two rows use L'(-a)~a, L''(-a)~-1, L'''(-a)~2/a^3.
  for (const auto& row :
       {std::array<double, 4>{-1e4, 10000.000099999998, -0.9999999900000006,
                              1.99999976000003e-12},
        {-1e8, 100000000.00000001, -0.9999999999999999, 1.9999999999999976e-24},
        {-1e10, 1e10, -1, 2e-30},
        {-1e100, 1e100, -1, 2e-300},
        {-1e308, 1e308, -1, 0}}) {
    SCOPED_TRACE(row[0]);
    check_tail_derivatives([](const auto& z) { return std_normal_lcdf(z); },
                           row[0], {row[1], row[2], row[3]});
    check_tail_derivatives(
        [](const auto& z) { return normal_lcdf(z, 0.0, 1.0); }, row[0],
        {row[1], row[2], row[3]});
    check_tail_derivatives([](const auto& z) { return std_normal_lccdf(z); },
                           -row[0], {-row[1], row[2], -row[3]});
    check_tail_derivatives(
        [](const auto& z) { return normal_lccdf(z, 0.0, 1.0); }, -row[0],
        {-row[1], row[2], -row[3]});
    check_tail_derivatives(
        [](const auto& z) {
          return internal::std_normal_lcdf_value_grad<false>(z).first;
        },
        row[0], {row[1], row[2], row[3]});
  }
}

TEST_F(AgradRev, std_normal_extreme_tail_vector_derivatives) {
  using namespace stan::math;
  const Eigen::Array<double, 5, 1> z(-1e8, -1e100, -1e308, 0, 1e308);
  const Eigen::Array<double, 5, 1> expected(2e-24, 2e-300, 0,
                                            0.21801361414499016, 0);
  Eigen::Array<fvar<fvar<var>>, 5, 1> x;
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    x(i).val_.val_ = z(i);
    x(i).val_.d_ = 1;
    x(i).d_.val_ = 1;
  }
  auto y = std_normal_lcdf(x);
  y.d_.d_.grad();
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    EXPECT_NEAR(expected(i), x(i).val_.val_.adj(),
                1e-12 * std::abs(expected(i)));
  }
}
