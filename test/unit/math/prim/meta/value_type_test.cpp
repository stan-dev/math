#include <test/unit/util.hpp>
#include <stan/math/prim/meta.hpp>
#include <unsupported/Eigen/KroneckerProduct>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, value_type_vector) {
  using stan::value_type;
  using std::vector;

  EXPECT_SAME_TYPE(vector<double>::value_type,
                   value_type<vector<double>>::type);

  EXPECT_SAME_TYPE(vector<double>::value_type,
                   value_type<const vector<double>>::type);

  EXPECT_SAME_TYPE(vector<vector<int>>::value_type,
                   value_type<vector<vector<int>>>::type);

  EXPECT_SAME_TYPE(vector<vector<int>>::value_type,
                   value_type<const vector<vector<int>>>::type);
}

TEST(MathMetaPrim, value_type_matrix) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::value_type;

  EXPECT_SAME_TYPE(Eigen::MatrixXd::Scalar, value_type<Eigen::MatrixXd>::type);

  EXPECT_SAME_TYPE(Eigen::VectorXd::Scalar, value_type<Eigen::VectorXd>::type);

  EXPECT_SAME_TYPE(Eigen::RowVectorXd::Scalar,
                   value_type<Eigen::RowVectorXd>::type);

  EXPECT_SAME_TYPE(Eigen::RowVectorXd,
                   value_type<std::vector<Eigen::RowVectorXd>>::type);
}

TEST(MathMetaPrim, value_type_kronecker) {
  Eigen::Matrix<double, 2, 2> A;
  const auto B
      = Eigen::kroneckerProduct(A, Eigen::Matrix<double, 2, 2>::Identity());
  Eigen::Matrix<double, 4, 1> C = Eigen::Matrix<double, 4, 1>::Random(4, 1);
  EXPECT_TRUE((std::is_same<double, stan::value_type_t<decltype(B)>>::value));
  Eigen::MatrixXd D = B * C;
}
