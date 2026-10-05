#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixProb_student_t_qf) {
  auto f = [](const auto& p, const auto& nu, const auto& sigma) {
    const auto& p_scaled = stan::math::inv_logit(p);
    const auto& nu_scaled = stan::math::exp(nu);
    const auto& sigma_scaled = stan::math::exp(sigma);
    return stan::math::sum(
        stan::math::student_t_qf(p_scaled, nu_scaled, 0, sigma_scaled));
  };

  stan::test::expect_ad(f, 0.3, 0.5, 3);
  stan::test::expect_ad(f, 0.8, 0.5, 0.1);
  stan::test::expect_ad(f, 0.1, 3, 3);

  Eigen::VectorXd p{{0.3, 0.8, 0.1}};

  Eigen::VectorXd nu{{0.5, 0.5, 3}};

  Eigen::VectorXd sigma{{3, 0.1, 3}};

  stan::test::expect_ad(f, p, nu, sigma);
}
