#ifndef TEST_UNIT_MATH_MIX_FUN_MULTIPLY_OPERATOR_UTIL_HPP
#define TEST_UNIT_MATH_MIX_FUN_MULTIPLY_OPERATOR_UTIL_HPP

#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <tuple>
#include <type_traits>

namespace stan {
namespace test {
namespace multiply_operator {

using stan::math::fvar;
using stan::math::var;

/**
 * The autodiff scalars `expect_ad` instantiates, for which Eigen products
 * must be evaluated by Stan's product hooks.
 */
using ad_scalars = std::tuple<var, fvar<double>, fvar<var>, fvar<fvar<double>>,
                              fvar<fvar<var>>>;

template <typename T, typename X>
struct promote {
  using type = stan::math::promote_scalar_t<T, X>;
};
template <typename T>
struct promote<T, double> {
  using type = T;
};
template <typename T>
struct promote<T, std::complex<double>> {
  using type = std::complex<T>;
};
template <typename T, typename S, int R, int C, int O, int MR, int MC>
struct promote<T, Eigen::Matrix<std::complex<S>, R, C, O, MR, MC>> {
  using type = Eigen::Matrix<std::complex<T>, R, C, O, MR, MC>;
};
template <typename T, typename X>
using promote_t = typename promote<T, X>::type;

/**
 * Expect Eigen to evaluate `xpr(a, b)` with a Stan product hook when either
 * or both arguments are promoted to autodiff scalar `T`. `xpr` returns the
 * unevaluated product expression; it is only used in `decltype`.
 */
template <typename T, typename F, typename A, typename B>
void expect_hooked_for(const F& xpr, const A&, const B&) {
  using a_t = promote_t<T, A>;
  using b_t = promote_t<T, B>;
  using both_t = std::decay_t<decltype(
      xpr(std::declval<const a_t&>(), std::declval<const b_t&>()))>;
  using left_t = std::decay_t<decltype(
      xpr(std::declval<const a_t&>(), std::declval<const B&>()))>;
  using right_t = std::decay_t<decltype(
      xpr(std::declval<const A&>(), std::declval<const b_t&>()))>;
  EXPECT_TRUE(stan::math::internal::is_stan_ad_product_hook_v<both_t>)
      << "both " << typeid(T).name();
  EXPECT_TRUE(stan::math::internal::is_stan_ad_product_hook_v<left_t>)
      << "left " << typeid(T).name();
  EXPECT_TRUE(stan::math::internal::is_stan_ad_product_hook_v<right_t>)
      << "right " << typeid(T).name();
}

/**
 * `expect_hooked_for` for every autodiff scalar, and expect the
 * all-`double` expression to stay with Eigen.
 */
template <typename F, typename A, typename B>
void expect_hooked(const F& xpr, const A& a, const B& b) {
  std::apply(
      [&](auto... t) { (expect_hooked_for<decltype(t)>(xpr, a, b), ...); },
      ad_scalars{});
  using double_t = std::decay_t<decltype(xpr(a, b))>;
  EXPECT_FALSE(stan::math::internal::is_stan_ad_product_hook_v<double_t>);
}

/**
 * Expect Eigen's assignment of a product to a destination with autodiff
 * scalar `T` (`=`, `+=`, `-=`, `noalias()`) to use Stan's hook, which
 * `Eigen::internal::Assignment` otherwise bypasses.
 */
template <typename T>
void expect_assignment_hooked_for() {
  using mat = Eigen::Matrix<T, -1, -1>;
  using product = Eigen::Product<mat, mat>;
  namespace ei = Eigen::internal;
  EXPECT_TRUE((std::is_base_of<stan::math::internal::ad_product_assignment<
                                   mat, mat, mat, 0, ei::assign_op<T, T>>,
                               ei::Assignment<mat, product, ei::assign_op<T, T>,
                                              ei::Dense2Dense>>::value));
  EXPECT_TRUE(
      (std::is_base_of<stan::math::internal::ad_product_assignment<
                           mat, mat, mat, 0, ei::add_assign_op<T, T>>,
                       ei::Assignment<mat, product, ei::add_assign_op<T, T>,
                                      ei::Dense2Dense>>::value));
  EXPECT_TRUE(
      (std::is_base_of<stan::math::internal::ad_product_assignment<
                           mat, mat, mat, 0, ei::sub_assign_op<T, T>>,
                       ei::Assignment<mat, product, ei::sub_assign_op<T, T>,
                                      ei::Dense2Dense>>::value));
}

inline void expect_assignment_hooked() {
  std::apply(
      [](auto... t) { (expect_assignment_hooked_for<decltype(t)>(), ...); },
      ad_scalars{});
}

inline Eigen::MatrixXd test_matrix(int rows, int cols, double offset) {
  Eigen::MatrixXd x(rows, cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      x(i, j) = 0.3 * (i + 1) - 0.2 * (j + 1) + offset;
    }
  }
  return x;
}

}  // namespace multiply_operator
}  // namespace test
}  // namespace stan

#endif
