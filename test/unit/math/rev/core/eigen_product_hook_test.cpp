#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/ad_tape_util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

// Eigen's `operator*` on `var` matrices is evaluated by Stan's free functions
// (multiply, elt_multiply, scale_multiply, multiply_add). Each test checks
// that an Eigen expression records the same nodes on the reverse-mode stack,
// and has the same value, as the direct call to the function it should use.

namespace {

using stan::math::var;
template <typename S, int R = Eigen::Dynamic, int C = Eigen::Dynamic>
using mat = Eigen::Matrix<S, R, C>;
template <typename S>
using vec = Eigen::Matrix<S, Eigen::Dynamic, 1>;
template <typename S>
using row_vec = Eigen::Matrix<S, 1, Eigen::Dynamic>;

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

/**
 * Expect `via_eigen()` to record the same chain and nochain varis as
 * `direct()` and to return the same values.
 */
template <typename Direct, typename ViaEigen>
void expect_same_as_direct(const Direct& direct, const ViaEigen& via_eigen) {
  const stan::test::tape_delta direct_delta;
  const auto expected = direct();
  const std::size_t direct_chain = direct_delta.chain();
  const std::size_t direct_nochain = direct_delta.nochain();
  const stan::test::tape_delta via_eigen_delta;
  const auto actual = via_eigen();
  EXPECT_EQ(via_eigen_delta.chain(), direct_chain);
  EXPECT_EQ(via_eigen_delta.nochain(), direct_nochain);
  EXPECT_MATRIX_NEAR(stan::math::value_of(actual),
                     stan::math::value_of(expected), 1e-12);
}

}  // namespace

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

TEST_F(AgradRev, eigen_product_hook_matrix_products_use_multiply) {
  mat<var> a = make_mat<var>(5, 4, 0.3);
  mat<var> b = make_mat<var>(4, 3, -0.2);
  const mat<double> b_d = make_mat<double>(4, 3, -0.2);
  vec<var> v = make_mat<var>(4, 1, 0.1);
  row_vec<var> rv = make_mat<var>(1, 5, 0.2);
  expect_same_as_direct([&] { return stan::math::multiply(a, b); },
                        [&] { return mat<var>(a * b); });
  expect_same_as_direct([&] { return stan::math::multiply(a, b_d); },
                        [&] { return mat<var>(a * b_d); });
  expect_same_as_direct(
      [&] { return stan::math::multiply(b_d.transpose(), a.transpose()); },
      [&] { return mat<var>(b_d.transpose() * a.transpose()); });
  expect_same_as_direct([&] { return stan::math::multiply(a, v); },
                        [&] { return vec<var>(a * v); });
  expect_same_as_direct([&] { return stan::math::multiply(v, rv); },
                        [&] { return mat<var>(v * rv); });
  expect_same_as_direct([&] { return stan::math::multiply(a, b); },
                        [&] { return mat<var>(a.lazyProduct(b)); });
  expect_same_as_direct(
      [&] {
        mat<var> c(5, 3);
        c.noalias() = stan::math::multiply(a, b);
        return c;
      },
      [&] {
        mat<var> c(5, 3);
        c.noalias() = a * b;
        return c;
      });
}

TEST_F(AgradRev, eigen_product_hook_fixed_size_products_use_multiply) {
  mat<var, 2, 2> f = make_mat<var>(2, 2, 0.4);
  const mat<double, 2, 2> g = make_mat<double>(2, 2, 0.1);
  mat<var, 4, 4> h = make_mat<var>(4, 4, -0.3);
  expect_same_as_direct([&] { return stan::math::multiply(f, g); },
                        [&] { return mat<var, 2, 2>(f * g); });
  expect_same_as_direct([&] { return stan::math::multiply(h, h); },
                        [&] { return mat<var, 4, 4>(h * h); });
}

TEST_F(AgradRev, eigen_product_hook_inner_product_uses_multiply) {
  row_vec<var> rv = make_mat<var>(1, 4, 0.2);
  vec<var> v = make_mat<var>(4, 1, 0.1);
  const stan::test::tape_delta direct_delta;
  const var expected = stan::math::multiply(rv, v);
  const std::size_t direct_total = direct_delta.total();
  const stan::test::tape_delta via_eigen_delta;
  const var actual = (rv * v)(0, 0);
  EXPECT_EQ(via_eigen_delta.total(), direct_total);
  EXPECT_NEAR(actual.val(), expected.val(), 1e-12);
}

TEST_F(AgradRev, eigen_product_hook_cwise_products_use_elt_multiply) {
  mat<var> x = make_mat<var>(4, 3, 0.3);
  mat<var> y = make_mat<var>(4, 3, -0.4);
  const mat<double> y_d = make_mat<double>(4, 3, -0.4);
  expect_same_as_direct([&] { return stan::math::elt_multiply(x, y); },
                        [&] { return mat<var>(x.array() * y.array()); });
  expect_same_as_direct([&] { return stan::math::elt_multiply(x, y_d); },
                        [&] { return mat<var>(x.cwiseProduct(y_d)); });
}

TEST_F(AgradRev, eigen_product_hook_scalar_products_use_multiply) {
  const mat<double> a = make_mat<double>(3, 2, 0.1);
  mat<var> a_v = make_mat<var>(3, 2, 0.1);
  var s = 1.5;
  expect_same_as_direct([&] { return stan::math::multiply(s, a); },
                        [&] { return mat<var>(s * a); });
  expect_same_as_direct([&] { return stan::math::multiply(s, a); },
                        [&] { return mat<var>(a * s); });
  expect_same_as_direct([&] { return stan::math::multiply(2.0, a_v); },
                        [&] { return mat<var>(a_v * 2); });
  expect_same_as_direct([&] { return stan::math::multiply(s, a_v); },
                        [&] { return mat<var>((a_v.array() * s).matrix()); });
}

TEST_F(AgradRev, eigen_product_hook_scaled_products_use_scale_multiply) {
  const mat<double> a = make_mat<double>(3, 4, 0.1);
  mat<var> b = make_mat<var>(4, 2, -0.3);
  var s = 0.7;
  expect_same_as_direct([&] { return stan::math::scale_multiply(s, a, b); },
                        [&] { return mat<var>(s * (a * b)); });
  expect_same_as_direct([&] { return stan::math::scale_multiply(s, a, b); },
                        [&] { return mat<var>((a * b) * s); });
  expect_same_as_direct([&] { return stan::math::scale_multiply(s, a, b); },
                        [&] { return mat<var>((s * a) * b); });
  expect_same_as_direct([&] { return stan::math::scale_multiply(s, a, b); },
                        [&] { return mat<var>(a * (b * s)); });
}

TEST_F(AgradRev, eigen_product_hook_assignments_use_multiply_add) {
  mat<var> a = make_mat<var>(3, 3, 0.3);
  mat<var> b = make_mat<var>(3, 3, -0.1);
  mat<var> c = make_mat<var>(3, 3, 0.05);
  expect_same_as_direct([&] { return stan::math::multiply_add(a, b, c); },
                        [&] {
                          mat<var> d = c;
                          d.noalias() += a * b;
                          return d;
                        });
  expect_same_as_direct([&] { return stan::math::multiply_add(a, b, c); },
                        [&] { return mat<var>(c + a * b); });
  expect_same_as_direct(
      [&] { return stan::math::subtract(c, stan::math::multiply(a, b)); },
      [&] {
        mat<var> d = c;
        d.noalias() -= a * b;
        return d;
      });
  // Without noalias() Eigen assumes the product may alias the destination,
  // so it evaluates the product into a temporary first and then adds it
  // coefficient by coefficient.
  expect_same_as_direct(
      [&] {
        mat<var> d = c;
        d += mat<var>(stan::math::multiply(a, b));
        return d;
      },
      [&] {
        mat<var> d = c;
        d += a * b;
        return d;
      });
}

TEST_F(AgradRev, eigen_product_hook_aliasing) {
  mat<var> a = make_mat<var>(3, 3, 0.3);
  mat<var> b = make_mat<var>(3, 3, -0.1);
  mat<var> e = a;
  e = e * b;
  EXPECT_MATRIX_NEAR(e.val(), a.val() * b.val(), 1e-12);
  mat<var> x = make_mat<var>(3, 1, 0.2);
  mat<var> x0 = x;
  x = x.cwiseProduct(x);
  EXPECT_MATRIX_NEAR(x.val(), x0.val().cwiseProduct(x0.val()), 1e-12);
}

TEST_F(AgradRev, eigen_product_hook_gradients) {
  const Eigen::Index n = 3, k = 4, m = 2;
  mat<var> a = make_mat<var>(n, k, 0.3);
  mat<var> b = make_mat<var>(k, m, -0.2);
  var s = 0.5;
  mat<var> c = s * (a * b);
  stan::math::sum(c).grad();
  const mat<double> ones = mat<double>::Ones(n, m);
  EXPECT_NEAR(s.adj(), (a.val() * b.val()).sum(), 1e-12);
  EXPECT_MATRIX_NEAR(a.adj(), 0.5 * ones * b.val().transpose(), 1e-12);
  EXPECT_MATRIX_NEAR(b.adj(), 0.5 * a.val().transpose() * ones, 1e-12);
}
