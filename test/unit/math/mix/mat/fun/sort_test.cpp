#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

void expect_sort(const std::vector<double>& sv,
                 const stan::test::ad_tolerances& tols
                 = stan::test::ad_tolerances()) {
  auto f_asc = [](const auto& x) { return stan::math::sort_asc(x); };
  auto f_desc = [](const auto& x) { return stan::math::sort_desc(x); };

  Eigen::VectorXd v(sv.size());
  Eigen::RowVectorXd rv(sv.size());
  for (size_t i = 0; i < sv.size(); ++i) {
    v(i) = sv[i];
    rv(i) = sv[i];
  }
  stan::test::expect_ad(tols, f_asc, sv);
  stan::test::expect_ad(tols, f_asc, v);
  stan::test::expect_ad(tols, f_asc, rv);
  stan::test::expect_ad(tols, f_desc, sv);
  stan::test::expect_ad(tols, f_desc, v);
  stan::test::expect_ad(tols, f_desc, rv);
}

TEST(MathMixMatFun, sort_asc_and_sort_desc) {
  std::vector<double> a;
  expect_sort(a);

  std::vector<double> b{1};
  expect_sort(b);

  std::vector<double> c1{1, 2};
  std::vector<double> c2{2, 1};
  expect_sort(c1);
  expect_sort(c2);

  std::vector<double> d1{1, 2, 3};
  std::vector<double> d2{1, 3, 2};
  std::vector<double> d3{2, 1, 3};
  std::vector<double> d4{2, 3, 1};
  std::vector<double> d5{3, 1, 2};
  std::vector<double> d6{3, 2, 1};
  expect_sort(d1);
  expect_sort(d2);
  expect_sort(d3);
  expect_sort(d4);
  expect_sort(d5);
  expect_sort(d6);

  // extreme tols to autodiff bigger and smaller numbers
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e0;
  tols.gradient_fvar_grad_ = 1e0;
  tols.hessian_grad_ = 1e0;
  tols.hessian_hessian_ = 1e5;
  tols.hessian_fvar_grad_ = 1e5;
  tols.hessian_fvar_hessian_ = 1e5;

  std::vector<double> w{-10.1, 2.12, 3.102};
  expect_sort(w, tols);

  std::vector<double> v{1.1e-6, -2.3, 31.1, 1, -10.1};
  expect_sort(v, tols);

  // lists with duplicates
  std::vector<double> e{1.1, 2.2, 33.1, -12.1, 33.1};
  expect_sort(e, tols);

  std::vector<double> u{1, -33.1, 2.1, -33.1};
  expect_sort(u, tols);

  // nan (exception)
  std::vector<double> h{1, std::numeric_limits<double>::quiet_NaN(), 3};
  expect_sort(h);

  // infinity
  std::vector<double> g{2, std::numeric_limits<double>::infinity(), 3, 4};
  expect_sort(g);
}
