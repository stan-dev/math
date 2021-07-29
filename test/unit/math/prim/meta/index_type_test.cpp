#include <stan/math/prim/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, index_type_vector) {
  using stan::math::index_type;
  using std::vector;

  EXPECT_SAME_TYPE(vector<double>::size_type,
                   index_type<vector<double> >::type);

  EXPECT_SAME_TYPE(vector<double>::size_type,
                   index_type<const vector<double> >::type);

  EXPECT_SAME_TYPE(vector<vector<int> >::size_type,
                   index_type<vector<vector<int> > >::type);

  EXPECT_SAME_TYPE(vector<vector<int> >::size_type,
                   index_type<const vector<vector<int> > >::type);
}

TEST(MathMetaPrim, index_type_matrices) {
  using stan::math::index_type;

  EXPECT_SAME_TYPE(Eigen::MatrixXd::Index, index_type<Eigen::MatrixXd>::type);

  EXPECT_SAME_TYPE(Eigen::VectorXd::Index, index_type<Eigen::VectorXd>::type);

  EXPECT_SAME_TYPE(Eigen::RowVectorXd::Index,
                   index_type<Eigen::RowVectorXd>::type);
}
