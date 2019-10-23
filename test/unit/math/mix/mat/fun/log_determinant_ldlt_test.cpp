#include <test/unit/math/test_ad.hpp>
// non-framework tests require following includes
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>

template <typename T, int R, int C>
stan::math::LDLT_factor<T, R, C> ldlt_factor(const Eigen::Matrix<T, R, C>& x) {
  stan::math::LDLT_factor<T, R, C> ldlt;
  ldlt.compute(x);
  return ldlt;
}

TEST(MathMixMatFun, logDeterminantLdlt) {
  auto f = [](const auto& x) {
    auto z = ((x + x.transpose()) * 0.5).eval();
    auto y = ldlt_factor(x);
    return stan::math::log_determinant_ldlt(y);
  };

  Eigen::MatrixXd a(2, 2);
  a << 3, 0, 0, 4;
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd c(2, 2);
  c << 1, 0, 0, 3;
  stan::test::expect_ad(f, c);

  for (const auto& rho : std::vector<double>{0, 0.9}) {
    for (int k = 1; k < 4; ++k) {
      Eigen::MatrixXd y(k, k);
      for (int i = 0; i < k; ++i) {
        y(i, i) = 0;
        for (int j = 0; j < i; ++j) {
          y(i, j) = std::pow(rho, std::fabs(i - j));
          y(j, i) = y(i, j);
        }
      }
      stan::test::expect_ad(f, y);
    }
  }

  // TODO(carpenter): figure out why this example requires loose tolerance
  // original derivative test for this case is included below
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1;  // default is tiny

  Eigen::MatrixXd b(2, 2);
  b << 2, 1, 1, 3;
  stan::test::expect_ad(tols, f, b);
}

TEST(MathMixMatFun, logDeterminantLdltSpecial) {
  using stan::math::determinant;
  using stan::math::fabs;
  using stan::math::log;
  using stan::math::matrix_v;

  // expected from auto-diff/Eigen
  AVEC x1 = createAVEC(2, 1, 1, 3);
  matrix_v v1(2, 2);
  v1 << x1[0], x1[1], x1[2], x1[3];
  AVAR det1 = log(fabs(v1.determinant()));
  std::vector<double> g1;
  det1.grad(x1, g1);

  stan::math::LDLT_factor<stan::math::var, -1, -1> ldlt_v;
  AVEC x2 = createAVEC(2, 1, 1, 3);
  matrix_v v2(2, 2);
  v2 << x2[0], x2[1], x2[2], x2[3];
  ldlt_v.compute(v2);
  ASSERT_TRUE(ldlt_v.success());
  AVAR det2 = log_determinant_ldlt(ldlt_v);
  std::vector<double> g2;
  det2.grad(x2, g2);

  EXPECT_FLOAT_EQ(det1.val(), det2.val());
  EXPECT_EQ(g1.size(), g2.size());
  for (size_t i = 0; i < g1.size(); ++i)
    EXPECT_FLOAT_EQ(g1[i], g2[i]);
}
