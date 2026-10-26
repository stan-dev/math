#include <test/unit/math/test_ad.hpp>
#include <vector>

namespace {

/**
 * Runs expect_ad at points in each region of grad_reg_inc_beta:
 *   (1, 1, 0.4), (2, 3, 0.25), (0.5, 0.5, 0.3): power series
 *   (10, 1, 0.5): power series, integer b
 *   (2, 3, 0.6): power series through the symmetry relation
 *   (20, 20, 0.3), (5, 30, 0.1): continued fraction
 *   (20, 20, 0.7), (5, 30, 0.3): continued fraction through the symmetry
 *     relation
 *   (10, 1, 0.7), (4, 4, 0.5): continued fraction, integer b, where the
 *     partial numerator d_{2b} is 0 but its b derivative is not
 */
template <typename F>
void expect_ad_grad_reg_inc_beta(const F& f) {
  stan::test::expect_ad(f, 1.0, 1.0, 0.4);
  stan::test::expect_ad(f, 2.0, 3.0, 0.25);
  stan::test::expect_ad(f, 0.5, 0.5, 0.3);
  stan::test::expect_ad(f, 10.0, 1.0, 0.5);
  stan::test::expect_ad(f, 2.0, 3.0, 0.6);
  stan::test::expect_ad(f, 20.0, 20.0, 0.3);
  stan::test::expect_ad(f, 5.0, 30.0, 0.1);
  stan::test::expect_ad(f, 20.0, 20.0, 0.7);
  stan::test::expect_ad(f, 5.0, 30.0, 0.3);
  stan::test::expect_ad(f, 10.0, 1.0, 0.7);
  stan::test::expect_ad(f, 4.0, 4.0, 0.5);
}

}  // namespace

TEST(mathMixScalFun, grad_reg_inc_beta_g1) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::digamma;
    return stan::math::grad_reg_inc_beta<stan::math::GRAD_RIB::G1>(
        a, b, z, digamma(a), 0.0, digamma(a + b), 0.0);
  };
  expect_ad_grad_reg_inc_beta(f);
}

TEST(mathMixScalFun, grad_reg_inc_beta_g2) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::digamma;
    return stan::math::grad_reg_inc_beta<stan::math::GRAD_RIB::G2>(
        a, b, z, 0.0, digamma(b), digamma(a + b), 0.0);
  };
  expect_ad_grad_reg_inc_beta(f);
}

TEST(mathMixScalFun, grad_reg_inc_beta_both) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::digamma;
    auto grads = stan::math::grad_reg_inc_beta<stan::math::GRAD_RIB::BOTH>(
        a, b, z, digamma(a), digamma(b), digamma(a + b), 0.0);
    Eigen::Matrix<decltype(grads.first), Eigen::Dynamic, 1> ret(2);
    ret << grads.first, grads.second;
    return ret;
  };
  expect_ad_grad_reg_inc_beta(f);
}

// the student_t_cdf call: a = nu / 2, b = 1 / 2 fixed, z = 1 - r
TEST(mathMixScalFun, grad_reg_inc_beta_student) {
  auto f = [](const auto& a, const auto& /* b */, const auto& z) {
    using stan::math::digamma;
    return stan::math::grad_reg_inc_beta<stan::math::GRAD_RIB::G1>(
        0.5 * a, 0.5, 1.0 - z, digamma(0.5 * a), 0.0, digamma(0.5 + 0.5 * a),
        0.0);
  };
  expect_ad_grad_reg_inc_beta(f);
}

/**
 * The overload writing both partials through out-parameters takes every
 * argument as the same type.
 */
TEST(mathMixScalFun, grad_reg_inc_beta_out_params) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::digamma;
    using T = stan::return_type_t<decltype(a), decltype(b), decltype(z)>;
    const T a_t = a;
    const T b_t = b;
    const T z_t = z;
    T g1 = 0;
    T g2 = 0;
    stan::math::grad_reg_inc_beta(g1, g2, a_t, b_t, z_t, digamma(a_t),
                                  digamma(b_t), digamma(a_t + b_t), T(0));
    Eigen::Matrix<T, Eigen::Dynamic, 1> ret(2);
    ret << g1, g2;
    return ret;
  };
  expect_ad_grad_reg_inc_beta(f);
}
