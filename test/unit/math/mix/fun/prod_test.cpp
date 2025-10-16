#include <test/unit/math/test_ad.hpp>
#include <vector>

template <typename T>
void expect_prod(const T& m) {
  auto f = [](const auto& x) { return stan::math::prod(x); };
  Eigen::VectorXd v(m.size());
  for (int i = 0; i < m.size(); ++i) {
    v(i) = m(i);
  }
  Eigen::RowVectorXd rv = v;
  stan::test::expect_ad(f, v);
  stan::test::expect_ad(f, rv);
  stan::test::expect_ad(f, m);
}

TEST(mathMixMatFun, prod) {
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m11(1, 1);
  m11 << 2;
  Eigen::MatrixXd m22(2, 2);
  m22 << 2, 3, 2, 3;
  Eigen::MatrixXd m23(2, 3);
  m23 << 10, 20, 30, -10, -30, -100;
  for (const auto& m : std::vector<Eigen::MatrixXd>{m00, m11, m22, m23}) {
    expect_prod(m);
  }
}
