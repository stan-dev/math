#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fma_vector) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    auto ret = stan::math::fma(x1, x2, x3).eval();
    using ret_t = std::decay_t<decltype(ret)>;
    constexpr bool is_correct_return_type
        = stan::is_var_matrix<ret_t>::value
          || stan::is_eigen_dense_base<ret_t>::value;
    static_assert(is_correct_return_type,
                  "Matrix type was supposed to be returned.");
    return ret;
  };

  double xd = 1.0;
  Eigen::VectorXd xv(2);
  xv << 1.0, 2.0;

  double yd = 2.0;
  Eigen::VectorXd yv(2);
  yv << 2.0, -3.0;

  double zd = 3.0;
  Eigen::VectorXd zv(2);
  zv << -1.0, 2.0;

  stan::test::expect_ad(f, xd, yd, zv);
  stan::test::expect_ad(f, xd, yv, zd);
  stan::test::expect_ad(f, xd, yv, zv);
  stan::test::expect_ad(f, xv, yd, zd);
  stan::test::expect_ad(f, xv, yd, zv);
  stan::test::expect_ad(f, xv, yv, zd);
  stan::test::expect_ad(f, xv, yv, zv);

  stan::test::expect_ad_matvar(f, xd, yd, zv);
  stan::test::expect_ad_matvar(f, xd, yv, zd);
  stan::test::expect_ad_matvar(f, xd, yv, zv);
  stan::test::expect_ad_matvar(f, xv, yd, zd);
  stan::test::expect_ad_matvar(f, xv, yd, zv);
  stan::test::expect_ad_matvar(f, xv, yv, zd);
  stan::test::expect_ad_matvar(f, xv, yv, zv);
}
