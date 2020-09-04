#include <stan/math/mix.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#define TEST_COMPARISON(SCALAR_TYPE, CONTAINER_TYPE1, CONTAINER_TYPE2, \
                        CONTAINER1, CONTAINER2, OP)                    \
  EXPECT_MATRIX_EQ(CONTAINER1 OP CONTAINER2,                           \
                   CONTAINER_TYPE1(CONTAINER1)                         \
                       OP CONTAINER_TYPE2(CONTAINER2));                \
  EXPECT_MATRIX_EQ(1 OP CONTAINER2,                                    \
                   SCALAR_TYPE(1) OP CONTAINER_TYPE2(CONTAINER2));     \
  EXPECT_MATRIX_EQ(CONTAINER1 OP 1,                                    \
                   CONTAINER_TYPE1(CONTAINER1) OP SCALAR_TYPE(1));

#define TEST_COMPARISON_COMBINATIONS(SCALAR_TYPE, CONTAINER_TYPE,      \
                                     CONTAINER_PLAIN_TYPE, CONTAINER1, \
                                     CONTAINER2, OP)                   \
  TEST_COMPARISON(double, CONTAINER_TYPE, CONTAINER_TYPE, CONTAINER1,  \
                  CONTAINER2, OP);                                     \
  TEST_COMPARISON(SCALAR_TYPE, CONTAINER_PLAIN_TYPE, CONTAINER_TYPE,   \
                  CONTAINER1, CONTAINER2, OP);                         \
  TEST_COMPARISON(SCALAR_TYPE, CONTAINER_TYPE, CONTAINER_PLAIN_TYPE,   \
                  CONTAINER1, CONTAINER2, OP);

#define TEST_COMPARISON_ALL_SHAPES(SCALAR_TYPE, OP)                            \
  {                                                                            \
    Eigen::ArrayXXd m1(3, 2);                                                  \
    m1 << -1, 2, 0.0, 0.5, 1, 1.5;                                             \
    Eigen::ArrayXXd m2(3, 2);                                                  \
    m2 << 2, 1, 0.0, -0.5, -1, 1.5;                                            \
    using T_m = Eigen::Array<SCALAR_TYPE, Eigen::Dynamic, Eigen::Dynamic>;     \
    TEST_COMPARISON_COMBINATIONS(SCALAR_TYPE, T_m, Eigen::ArrayXXd, m1, m2,    \
                                 OP);                                          \
    Eigen::ArrayXd v1(6);                                                   \
    v1 << -1, 2, 0.0, 0.5, 1, 1.5;                                             \
    Eigen::ArrayXd v2(6);                                                   \
    v2 << 2, 1, 0.0, -0.5, -1, 1.5;                                            \
    using T_v = Eigen::Array<SCALAR_TYPE, Eigen::Dynamic, 1>;                  \
    TEST_COMPARISON_COMBINATIONS(SCALAR_TYPE, T_v, Eigen::ArrayXd, v1, v2,     \
                                 OP);                                          \
    using T_rv_plain = Eigen::Array<double, 1, Eigen::Dynamic>;                \
    T_rv_plain rv1(6);                                                      \
    rv1 << -1, 2, 0.0, 0.5, 1, 1.5;                                            \
    T_rv_plain rv2(6);                                                      \
    rv2 << 2, 1, 0.0, -0.5, -1, 1.5;                                           \
    using T_rv = Eigen::Array<SCALAR_TYPE, 1, Eigen::Dynamic>;                 \
    TEST_COMPARISON_COMBINATIONS(SCALAR_TYPE, T_rv, T_rv_plain, rv1, rv2, OP); \
  }

#define TEST_ALL_COMPARISONS(NAME, SCALAR_TYPE)  \
  TEST(mixFun, eigen_comparison_##NAME) {        \
    TEST_COMPARISON_ALL_SHAPES(SCALAR_TYPE, <);  \
    TEST_COMPARISON_ALL_SHAPES(SCALAR_TYPE, <=); \
    TEST_COMPARISON_ALL_SHAPES(SCALAR_TYPE, >);  \
    TEST_COMPARISON_ALL_SHAPES(SCALAR_TYPE, >=); \
    TEST_COMPARISON_ALL_SHAPES(SCALAR_TYPE, ==); \
    TEST_COMPARISON_ALL_SHAPES(SCALAR_TYPE, !=); \
  }

TEST_ALL_COMPARISONS(var, stan::math::var)
TEST_ALL_COMPARISONS(fvar_double, stan::math::fvar<double>)
TEST_ALL_COMPARISONS(fvar_var, stan::math::fvar<stan::math::var>)
TEST_ALL_COMPARISONS(fvar_fvar_double,
                     stan::math::fvar<stan::math::fvar<double>>)
