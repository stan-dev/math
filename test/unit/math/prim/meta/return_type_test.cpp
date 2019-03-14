
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>






using stan::return_type;

TEST(MetaTraits, ReturnTypeDouble) {
  test::expect_same_type<double, return_type<double>::type>();
}

TEST(MetaTraits, ReturnTypeFloat) {
  test::expect_same_type<double, return_type<float>::type>();
}

TEST(MetaTraits, ReturnTypeInt) {
  test::expect_same_type<double, return_type<int>::type>();
}

TEST(MetaTraits, ReturnTypeScalarTenParams) {
  test::expect_same_type<double,
                         return_type<double, int, double, float, float, double,
                                     float, int, double, double>::type>();
}





using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using stan::return_type;
using std::vector;

TEST(MetaTraits_mat, ReturnTypeMatrixXd) {
  test::expect_same_type<double, return_type<MatrixXd>::type>();
}

TEST(MetaTraits_mat, ReturnTypeVectorXd) {
  test::expect_same_type<double, return_type<VectorXd>::type>();
}

TEST(MetaTraits_mat, ReturnTypeRowVectorXd) {
  test::expect_same_type<double, return_type<RowVectorXd>::type>();
}

TEST(MetaTraits_mat, ReturnTypeArray) {
  test::expect_same_type<double, return_type<vector<int> >::type>();
}
