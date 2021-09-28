#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <vector>

// test sum of first n numbers for sum of a
template <typename T>
void test_sum(stan::math::accumulator<T>& a, int n) {
  EXPECT_TRUE((n * (n + 1)) / 2 == a.sum());
}

TEST(AgradMixMatrixAccumulate, fvar_var) {
  using stan::math::accumulator;
  using stan::math::fvar;
  using stan::math::var;

  accumulator<fvar<var>> a;
  test_sum(a, 0);

  a.add(1.0);
  test_sum(a, 1);

  for (int i = 2; i <= 1000; ++i)
    a.add(i);
  test_sum(a, 1000);
}

TEST(AgradMixMatrixAccumulate, collection_fvar_var) {
  using stan::math::accumulator;
  using stan::math::fvar;
  using stan::math::matrix_fv;
  using stan::math::var;
  using stan::math::vector_fv;
  using std::vector;

  accumulator<fvar<var>> a;

  int pos = 0;
  test_sum(a, 0);

  vector<fvar<var>> v(10);
  for (size_t i = 0; i < 10; ++i)
    v[i] = pos++;
  a.add(v);
  test_sum(a, pos - 1);

  a.add(pos++);
  test_sum(a, pos - 1);

  double x = pos++;
  a.add(x);
  test_sum(a, pos - 1);

  vector<vector<fvar<var>>> ww(10);
  for (size_t i = 0; i < 10; ++i) {
    vector<fvar<var>> w(5);
    for (size_t n = 0; n < 5; ++n)
      w[n] = pos++;
    ww[i] = w;
  }
  a.add(ww);
  test_sum(a, pos - 1);

  matrix_fv m(5, 6);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 6; ++j)
      m(i, j) = pos++;
  a.add(m);
  test_sum(a, pos - 1);

  vector_fv mv(7);
  for (int i = 0; i < 7; ++i)
    mv(i) = pos++;
  a.add(mv);
  test_sum(a, pos - 1);

  vector<vector_fv> vvx(8);
  for (size_t i = 0; i < 8; ++i) {
    vector_fv vx(3);
    for (int j = 0; j < 3; ++j)
      vx(j) = pos++;
    vvx[i] = vx;
  }
  a.add(vvx);
  test_sum(a, pos - 1);
}

TEST(AgradMixMatrixAccumulate, fvar_fvar_var) {
  using stan::math::accumulator;
  using stan::math::fvar;
  using stan::math::var;

  accumulator<fvar<fvar<var>>> a;
  test_sum(a, 0);

  a.add(1.0);
  test_sum(a, 1);

  for (int i = 2; i <= 1000; ++i)
    a.add(i);
  test_sum(a, 1000);
}

TEST(AgradMixMatrixAccumulate, collection_fvar_fvar_var) {
  using stan::math::accumulator;
  using stan::math::fvar;
  using stan::math::matrix_ffv;
  using stan::math::var;
  using stan::math::vector_ffv;
  using std::vector;

  accumulator<fvar<fvar<var>>> a;

  int pos = 0;
  test_sum(a, 0);

  vector<fvar<fvar<var>>> v(10);
  for (size_t i = 0; i < 10; ++i)
    v[i] = pos++;
  a.add(v);
  test_sum(a, pos - 1);

  a.add(pos++);
  test_sum(a, pos - 1);

  int x = pos++;
  a.add(x);
  test_sum(a, pos - 1);

  vector<vector<fvar<fvar<var>>>> ww(10);
  for (size_t i = 0; i < 10; ++i) {
    vector<fvar<fvar<var>>> w(5);
    for (size_t n = 0; n < 5; ++n)
      w[n] = pos++;
    ww[i] = w;
  }
  a.add(ww);
  test_sum(a, pos - 1);

  matrix_ffv m(5, 6);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 6; ++j)
      m(i, j) = pos++;
  a.add(m);
  test_sum(a, pos - 1);

  vector_ffv mv(7);
  for (int i = 0; i < 7; ++i)
    mv(i) = pos++;
  a.add(mv);
  test_sum(a, pos - 1);

  vector<vector_ffv> vvx(8);
  for (size_t i = 0; i < 8; ++i) {
    vector_ffv vx(3);
    for (int j = 0; j < 3; ++j)
      vx(j) = pos++;
    vvx[i] = vx;
  }
  a.add(vvx);
  test_sum(a, pos - 1);
}

TEST(AgradMixMatrixAccumulate, var_matrix) {
  auto f = [](const auto& x) {
    using x_t = decltype(x);
    stan::math::accumulator<stan::scalar_type_t<x_t>> acc;
    acc.add(x);
    return acc.sum();
  };
  Eigen::MatrixXd x(2, 2);
  x << 1, -1, 10, 0;
  stan::test::expect_ad(f, x);
  stan::test::expect_ad_matvar(f, x);
  std::vector<Eigen::MatrixXd> x_vec;
  x_vec.push_back(x);
  x_vec.push_back(x);
  x_vec.push_back(x);
  stan::test::expect_ad(f, x_vec);
  stan::test::expect_ad_matvar(f, x_vec);
}
