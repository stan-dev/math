#include <stan/math/prim/meta.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/functor/reduce_sum_util.hpp>
#include <test/unit/math/rev/util.hpp>

#include <limits>
#include <vector>

// Reduce sum tests are broken up into four files to avoid windows compiler
// error

TEST_F(AgradRev, reduce_sum_eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(2);
  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST_F(AgradRev, reduce_sum_eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(2);

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST_F(AgradRev, reduce_sum_eigen_matrix_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(2, 2);

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST_F(AgradRev, reduce_sum_std_vector_std_vector_double_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<double>> arg(2, std::vector<double>(2, 10.0));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST_F(AgradRev, reduce_sum_std_vector_eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(2));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST_F(AgradRev, reduce_sum_std_vector_eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::RowVectorXd> arg(2, Eigen::RowVectorXd::Ones(2));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}
