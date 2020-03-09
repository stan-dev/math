#include <test/unit/math/test_ad.hpp>

template <typename T>
auto g1(const T& x) {
  auto x_cons = stan::math::simplex_constrain(x);
  auto x_free = stan::math::simplex_free(x_cons);
  return x_free;
}
template <typename T>
auto g2(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  auto x_cons = stan::math::simplex_constrain(x, lp);
  auto x_free = stan::math::simplex_free(x_cons);
  return x_free;
}
template <typename T>
auto g3(const T& x) {
  stan::scalar_type_t<T> lp = 0;
  auto x_cons = stan::math::simplex_constrain(x, lp);
  auto x_free = stan::math::simplex_free(x_cons);
  return lp;
}

template <typename T>
void expect_simplex_transform(const T& x) {
  auto f1 = [](const auto& x) { return g1(x); };
  auto f2 = [](const auto& x) { return g2(x); };
  auto f3 = [](const auto& x) { return g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad(f3, x);
}

TEST(MathMixMatFun, simplexTransformMat) {
  Eigen::VectorXd v0(0);
  expect_simplex_transform(v0);

  Eigen::VectorXd v1(1);
  v1 << 1;
  expect_simplex_transform(v1);

  Eigen::VectorXd v2(2);
  v2 << 3, -1;
  expect_simplex_transform(v2);

  Eigen::VectorXd v3(3);
  v3 << 2, 3, -1;
  expect_simplex_transform(v3);

  Eigen::VectorXd v4(4);
  v4 << 2, -1, 0, -1.1;
  expect_simplex_transform(v4);

  Eigen::VectorXd v5(5);
  v5 << 1, -3, 2, 0, -1;
  expect_simplex_transform(v5);
}

TEST(MathMixMatFun, simplexTransformVec) {
  std::vector<double> v0(0);
  expect_simplex_transform(v0);

  std::vector<double> v1({1});
  expect_simplex_transform(v1);

  std::vector<double> v2({3, -1});
  expect_simplex_transform(v2);

  std::vector<double> v3({-12, 3, -1.9});
  expect_simplex_transform(v3);

  std::vector<double> v4({-1, 0, -1.1, 0.5});
  expect_simplex_transform(v4);

  std::vector<double> v5({1, -3, 2, 0, -1});
  expect_simplex_transform(v5);
}
