#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

using stan::return_type;

TEST(MathMetaPrim, ReturnType_scalar) {
  test::expect_same_type<double, return_type<double>::type>();
  test::expect_same_type<double, return_type<float>::type>();
  test::expect_same_type<double, return_type<int>::type>();
  test::expect_same_type<double,
                         return_type<double, int, double, float, float, double,
                                     float, int, double, double>::type>();
}

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using stan::return_type;
using std::vector;

TEST(MathMetaPrim, ReturnType_non_scalar) {
  test::expect_same_type<double, return_type<MatrixXd>::type>();
  test::expect_same_type<double, return_type<VectorXd>::type>();
  test::expect_same_type<double, return_type<RowVectorXd>::type>();
  test::expect_same_type<double, return_type<vector<int> >::type>();
}
