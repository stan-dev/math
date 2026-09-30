#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/ad_tape_util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <complex>
#include <type_traits>

namespace {

using stan::math::var;
using stan::math::internal::is_stan_ad_product_hook_v;
using cvar = std::complex<var>;
template <typename S, int R = Eigen::Dynamic, int C = Eigen::Dynamic>
using mat = Eigen::Matrix<S, R, C>;
template <typename S>
using vec = Eigen::Matrix<S, Eigen::Dynamic, 1>;
template <typename S>
using row_vec = Eigen::Matrix<S, 1, Eigen::Dynamic>;
template <typename X, typename Y>
using product_t = std::decay_t<decltype(std::declval<const X&>()
                                        * std::declval<const Y&>())>;

template <typename S>
mat<S> make_mat(Eigen::Index rows, Eigen::Index cols, double offset) {
  mat<S> x(rows, cols);
  for (Eigen::Index j = 0; j < cols; ++j) {
    for (Eigen::Index i = 0; i < rows; ++i) {
      x(i, j) = 0.1 * (i + 1) - 0.07 * (j + 1) + offset;
    }
  }
  return x;
}

template <typename Derived>
var sum_of(const Eigen::MatrixBase<Derived>& x) {
  var total = 0;
  for (Eigen::Index j = 0; j < x.cols(); ++j) {
    for (Eigen::Index i = 0; i < x.rows(); ++i) {
      total += x.derived().coeff(i, j);
    }
  }
  return total;
}

}  // namespace

TEST_F(AgradRev, eigen_product_hook_selects_var_products) {
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<mat<var>, mat<double>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<mat<double>, mat<var>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<mat<var>, mat<var>>>));
  EXPECT_TRUE(
      (is_stan_ad_product_hook_v<product_t<mat<var, 2, 2>, mat<var, 2, 2>>>));
  EXPECT_TRUE((
      is_stan_ad_product_hook_v<product_t<mat<var, 4, 4>, mat<double, 4, 4>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<mat<var>, vec<double>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<row_vec<var>, vec<var>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<vec<var>, row_vec<var>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<mat<cvar>, mat<var>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<
               product_t<mat<std::complex<double>>, mat<cvar>>>));
  EXPECT_TRUE(
      (is_stan_ad_product_hook_v<
          std::decay_t<decltype(std::declval<const mat<var>&>().lazyProduct(
              std::declval<const mat<double>&>()))>>));
  EXPECT_FALSE(
      (is_stan_ad_product_hook_v<product_t<mat<double>, mat<double>>>));
  EXPECT_FALSE((is_stan_ad_product_hook_v<
                product_t<mat<double, 2, 2>, mat<double, 2, 2>>>));
}

TEST_F(AgradRev, eigen_product_hook_selects_var_cwise_products) {
  using arr = Eigen::Array<var, Eigen::Dynamic, Eigen::Dynamic>;
  using arr_d = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<arr, arr_d>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<arr_d, arr>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<arr, arr>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<var, mat<double>>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<mat<double>, var>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<mat<var>, int>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<product_t<double, mat<var>>>));
  EXPECT_TRUE(
      (is_stan_ad_product_hook_v<
          std::decay_t<decltype(std::declval<const mat<var>&>().cwiseProduct(
              std::declval<const mat<double>&>()))>>));
  EXPECT_FALSE((is_stan_ad_product_hook_v<product_t<arr_d, arr_d>>));
  EXPECT_FALSE((is_stan_ad_product_hook_v<product_t<double, mat<double>>>));
}

TEST_F(AgradRev, eigen_product_hook_double_products_unchanged) {
  const mat<double> a = make_mat<double>(3, 4, 0.2);
  const mat<double> b = make_mat<double>(4, 2, -0.1);
  const mat<double> c = a * b;
  const mat<double> c_lazy = a.lazyProduct(b);
  for (Eigen::Index j = 0; j < 2; ++j) {
    for (Eigen::Index i = 0; i < 3; ++i) {
      double expected = 0;
      for (Eigen::Index k = 0; k < 4; ++k) {
        expected += a(i, k) * b(k, j);
      }
      EXPECT_DOUBLE_EQ(c(i, j), expected);
      EXPECT_DOUBLE_EQ(c_lazy(i, j), expected);
    }
  }
  EXPECT_EQ(stan::test::tape_size(), 0);
}

TEST_F(AgradRev, eigen_product_hook_matmul_one_callback) {
  const Eigen::Index n = 5, k = 4, m = 3;
  mat<var> a = make_mat<var>(n, k, 0.3);
  const mat<double> b = make_mat<double>(k, m, -0.2);
  stan::test::tape_delta delta;
  mat<var> c = a * b;
  EXPECT_EQ(delta.chain(), 1);
  EXPECT_EQ(delta.nochain(), n * m);
  const mat<double> expected = a.val() * b;
  EXPECT_MATRIX_NEAR(c.val(), expected, 1e-12);
  sum_of(c).grad();
  // d sum(A * B) / dA = ones(n, m) * B^T
  const mat<double> grad_a = mat<double>::Ones(n, m) * b.transpose();
  EXPECT_MATRIX_NEAR(a.adj(), grad_a, 1e-12);
}

TEST_F(AgradRev, eigen_product_hook_matmul_var_var) {
  const Eigen::Index n = 3, k = 4, m = 2;
  mat<var> a = make_mat<var>(n, k, 0.3);
  mat<var> b = make_mat<var>(k, m, -0.2);
  stan::test::tape_delta delta;
  mat<var> c = a * b;
  EXPECT_EQ(delta.chain(), 1);
  EXPECT_EQ(delta.nochain(), n * m);
  sum_of(c).grad();
  EXPECT_MATRIX_NEAR(a.adj(), mat<double>::Ones(n, m) * b.val().transpose(),
                     1e-12);
  EXPECT_MATRIX_NEAR(b.adj(), a.val().transpose() * mat<double>::Ones(n, m),
                     1e-12);
}

TEST_F(AgradRev, eigen_product_hook_matmul_shapes) {
  vec<var> v = make_mat<var>(4, 1, 0.1);
  row_vec<var> rv = make_mat<var>(1, 4, 0.2);
  const mat<double> a = make_mat<double>(3, 4, 0.3);
  {
    stan::test::tape_delta delta;
    vec<var> av = a * v;
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), 3);
    EXPECT_MATRIX_NEAR(av.val(), a * v.val(), 1e-12);
  }
  {
    stan::test::tape_delta delta;
    var dot = (rv * v)(0, 0);
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), 1);
    EXPECT_NEAR(dot.val(), rv.val().dot(v.val()), 1e-12);
  }
  {
    stan::test::tape_delta delta;
    mat<var> outer = v * rv;
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), 16);
    EXPECT_MATRIX_NEAR(outer.val(), v.val() * rv.val(), 1e-12);
  }
  {
    mat<var, 2, 2> f = make_mat<var>(2, 2, 0.4);
    mat<double, 2, 2> g = make_mat<double>(2, 2, 0.1);
    stan::test::tape_delta delta;
    mat<var, 2, 2> fg = f * g;
    mat<var, 2, 2> fg_lazy = f.lazyProduct(g);
    EXPECT_EQ(delta.chain(), 2);
    EXPECT_EQ(delta.nochain(), 8);
    EXPECT_MATRIX_NEAR(fg.val(), f.val() * g, 1e-12);
    EXPECT_MATRIX_NEAR(fg_lazy.val(), f.val() * g, 1e-12);
  }
}

TEST_F(AgradRev, eigen_product_hook_assignment_forms) {
  mat<var> a = make_mat<var>(3, 3, 0.3);
  mat<var> b = make_mat<var>(3, 3, -0.1);
  mat<var> c = make_mat<var>(3, 3, 0.05);
  const mat<double> c0 = c.val();
  {
    stan::test::tape_delta delta;
    mat<var> d(3, 3);
    d.noalias() = a * b;
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), 9);
  }
  {
    stan::test::tape_delta delta;
    c += a * b;
    // one callback for the product plus one vari per coefficient for +=
    EXPECT_EQ(delta.chain(), 1 + 9);
    EXPECT_EQ(delta.nochain(), 9);
    EXPECT_MATRIX_NEAR(c.val(), c0 + a.val() * b.val(), 1e-12);
  }
  {
    // aliasing: the product is evaluated before the destination is written
    mat<var> e = a;
    e = e * b;
    EXPECT_MATRIX_NEAR(e.val(), a.val() * b.val(), 1e-12);
  }
}

TEST_F(AgradRev, eigen_product_hook_cwise_one_callback) {
  const Eigen::Index n = 4;
  vec<var> x = make_mat<var>(n, 1, 0.3);
  vec<var> y = make_mat<var>(n, 1, -0.4);
  stan::test::tape_delta delta;
  vec<var> z = (x.array() * y.array()).matrix();
  EXPECT_EQ(delta.chain(), 1);
  EXPECT_EQ(delta.nochain(), n);
  sum_of(z).grad();
  EXPECT_MATRIX_NEAR(x.adj(), y.val(), 1e-12);
  EXPECT_MATRIX_NEAR(y.adj(), x.val(), 1e-12);
}

TEST_F(AgradRev, eigen_product_hook_scalar_times_matrix) {
  const mat<double> a = make_mat<double>(3, 2, 0.1);
  var s = 1.5;
  {
    stan::test::tape_delta delta;
    mat<var> c = s * a;
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), 6);
    EXPECT_MATRIX_NEAR(c.val(), 1.5 * a, 1e-12);
    sum_of(c).grad();
    EXPECT_NEAR(s.adj(), a.sum(), 1e-12);
  }
  stan::math::set_zero_all_adjoints();
  mat<var> av = a;
  {
    stan::test::tape_delta delta;
    mat<var> c = av * 2;
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), 6);
    sum_of(c).grad();
    EXPECT_MATRIX_NEAR(av.adj(), mat<double>::Constant(3, 2, 2.0), 1e-12);
  }
}

TEST_F(AgradRev, eigen_product_hook_scaled_products_fuse) {
  const Eigen::Index n = 3, k = 4, m = 2;
  const mat<double> a = make_mat<double>(n, k, 0.1);
  mat<var> b = make_mat<var>(k, m, -0.3);
  var s = 0.7;
  const mat<double> expected = 0.7 * (a * b.val());
  {
    stan::test::tape_delta delta;
    mat<var> c = s * (a * b);
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), n * m);
    EXPECT_MATRIX_NEAR(c.val(), expected, 1e-12);
    sum_of(c).grad();
    EXPECT_NEAR(s.adj(), (a * b.val()).sum(), 1e-12);
    EXPECT_MATRIX_NEAR(b.adj(), 0.7 * a.transpose() * mat<double>::Ones(n, m),
                       1e-12);
  }
  {
    stan::test::tape_delta delta;
    mat<var> c = (a * b) * s;
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), n * m);
    EXPECT_MATRIX_NEAR(c.val(), expected, 1e-12);
  }
  {
    stan::test::tape_delta delta;
    mat<var> c = (s * a) * b;
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), n * m);
    EXPECT_MATRIX_NEAR(c.val(), expected, 1e-12);
  }
  {
    stan::test::tape_delta delta;
    mat<var> c = a * (b * s);
    EXPECT_EQ(delta.chain(), 1);
    EXPECT_EQ(delta.nochain(), n * m);
    EXPECT_MATRIX_NEAR(c.val(), expected, 1e-12);
  }
  {
    mat<var> d = make_mat<var>(n, m, 0.2);
    stan::test::tape_delta delta;
    mat<var> c = d + s * (a * b);
    // one product callback plus one vari per coefficient for the sum
    EXPECT_EQ(delta.chain(), 1 + n * m);
    EXPECT_EQ(delta.nochain(), n * m);
    EXPECT_MATRIX_NEAR(c.val(), d.val() + expected, 1e-12);
  }
}

TEST_F(AgradRev, eigen_product_hook_complex_matmul) {
  mat<cvar> a(2, 2);
  a << cvar(1.0, 0.5), cvar(-0.3, 0.2), cvar(0.4, -1.0), cvar(2.0, 0.1);
  mat<std::complex<double>> b(2, 1);
  b << std::complex<double>(0.5, -0.25), std::complex<double>(1.5, 0.75);
  stan::test::tape_delta delta;
  mat<cvar> c = a * b;
  EXPECT_EQ(delta.chain(), 1);
  EXPECT_EQ(delta.nochain(), 4);
  mat<std::complex<double>> a_val(2, 2);
  for (Eigen::Index j = 0; j < 2; ++j) {
    for (Eigen::Index i = 0; i < 2; ++i) {
      a_val(i, j) = {a(i, j).real().val(), a(i, j).imag().val()};
    }
  }
  const mat<std::complex<double>> expected = a_val * b;
  for (Eigen::Index i = 0; i < 2; ++i) {
    EXPECT_NEAR(c(i, 0).real().val(), expected(i, 0).real(), 1e-12);
    EXPECT_NEAR(c(i, 0).imag().val(), expected(i, 0).imag(), 1e-12);
  }
  // d re(c_0) / d re(a_00) = re(b_0), d re(c_0) / d im(a_00) = -im(b_0)
  c(0, 0).real().grad();
  EXPECT_NEAR(a(0, 0).real().adj(), b(0).real(), 1e-12);
  EXPECT_NEAR(a(0, 0).imag().adj(), -b(0).imag(), 1e-12);
  EXPECT_NEAR(a(1, 0).real().adj(), 0.0, 1e-12);
}
