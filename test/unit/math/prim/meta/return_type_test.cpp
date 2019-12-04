#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

using stan::return_type;

TEST(MetaTraitsPrimScal, ReturnTypeDouble) {
  test::expect_same_type<double, return_type<double>::type>();
}

TEST(MetaTraitsPrimScal, ReturnTypeFloat) {
  test::expect_same_type<double, return_type<float>::type>();
}

TEST(MetaTraitsPrimScal, ReturnTypeInt) {
  test::expect_same_type<double, return_type<int>::type>();
}

TEST(MetaTraitsPrimScal, ReturnTypeScalarTenParams) {
  test::expect_same_type<double,
                         return_type<double, int, double, float, float, double,
                                     float, int, double, double>::type>();
}

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using stan::return_type;
using std::vector;

TEST(MetaTraitsPrimMat, ReturnTypeMatrixXd) {
  test::expect_same_type<double, return_type<MatrixXd>::type>();
}

TEST(MetaTraitsPrimMat, ReturnTypeVectorXd) {
  test::expect_same_type<double, return_type<VectorXd>::type>();
}

TEST(MetaTraitsPrimMat, ReturnTypeRowVectorXd) {
  test::expect_same_type<double, return_type<RowVectorXd>::type>();
}

TEST(MetaTraitsPrimMat, ReturnTypeArray) {
  test::expect_same_type<double, return_type<vector<int> >::type>();
}
