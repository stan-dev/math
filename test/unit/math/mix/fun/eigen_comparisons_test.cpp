#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

template <typename Scalar, typename ContainerT1, typename ContainerT2,
          typename Op, typename Container1Plain, typename Container2Plain>
void test_comparison(Op operation, const Container1Plain& container1_plain,
                     const Container2Plain& container2_plain) {
  EXPECT_MATRIX_EQ(
      operation(container1_plain, container2_plain),
      operation(ContainerT1(container1_plain), ContainerT2(container2_plain)));
  EXPECT_MATRIX_EQ(operation(1, container2_plain),
                   operation(Scalar(1), ContainerT2(container2_plain)));
  EXPECT_MATRIX_EQ(operation(container1_plain, 1),
                   operation(ContainerT1(container1_plain), Scalar(1)));
}

template <typename Scalar, typename ContainerT, typename Op,
          typename ContainerPlain>
void test_comparison_combinations(Op operation,
                                  const ContainerPlain& container1_plain,
                                  const ContainerPlain& container2_plain) {
  test_comparison<double, ContainerT, ContainerT>(operation, container1_plain,
                                                  container2_plain);
  test_comparison<Scalar, ContainerPlain, ContainerT>(
      operation, container1_plain, container2_plain);
  test_comparison<Scalar, ContainerT, ContainerPlain>(
      operation, container1_plain, container2_plain);
}

template <typename Scalar, typename Op>
void test_comparison_all_shapes(Op operation) {
  Eigen::ArrayXXd m1(3, 2);
  m1 << -1, 2, 0.0, 0.5, 1, 1.5;
  Eigen::ArrayXXd m2(3, 2);
  m2 << 2, 1, 0.0, -0.5, -1, 1.5;
  using T_m = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  test_comparison_combinations<Scalar, T_m>(operation, m1, m2);

  Eigen::ArrayXd v1(6);
  v1 << -1, 2, 0.0, 0.5, 1, 1.5;
  Eigen::ArrayXd v2(6);
  v2 << 2, 1, 0.0, -0.5, -1, 1.5;
  using T_v = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
  test_comparison_combinations<Scalar, T_v>(operation, v1, v2);

  using T_rv_plain = Eigen::Array<double, 1, Eigen::Dynamic>;
  T_rv_plain rv1(6);
  rv1 << -1, 2, 0.0, 0.5, 1, 1.5;
  T_rv_plain rv2(6);
  rv2 << 2, 1, 0.0, -0.5, -1, 1.5;
  using T_rv = Eigen::Array<Scalar, 1, Eigen::Dynamic>;
  test_comparison_combinations<Scalar, T_rv>(operation, rv1, rv2);
}

template <typename Scalar>
void test_all_comparisons() {
  test_comparison_all_shapes<Scalar>(
      [](const auto& a, const auto& b) { return a < b; });
  test_comparison_all_shapes<Scalar>(
      [](const auto& a, const auto& b) { return a <= b; });
  test_comparison_all_shapes<Scalar>(
      [](const auto& a, const auto& b) { return a > b; });
  test_comparison_all_shapes<Scalar>(
      [](const auto& a, const auto& b) { return a >= b; });
  test_comparison_all_shapes<Scalar>(
      [](const auto& a, const auto& b) { return a == b; });
  test_comparison_all_shapes<Scalar>(
      [](const auto& a, const auto& b) { return a != b; });
}

TEST(mixFun, eigen_comparisons_var) { test_all_comparisons<stan::math::var>(); }

TEST(mixFun, eigen_comparisons_fvar) {
  test_all_comparisons<stan::math::fvar<double>>();
}

TEST(mixFun, eigen_comparisons_fvar_var) {
  test_all_comparisons<stan::math::fvar<stan::math::var>>();
}

TEST(mixFun, eigen_comparisons_fvar_fvar) {
  test_all_comparisons<stan::math::fvar<stan::math::fvar<double>>>();
}
