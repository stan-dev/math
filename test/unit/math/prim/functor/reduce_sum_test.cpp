#include <stan/math.hpp>
#include <test/unit/math/prim/functor/reduce_sum_util.hpp>
#include <gtest/gtest.h>

#include <vector>
#include <set>

TEST(StanMathPrim_reduce_sum, value) {
  using stan::math::test::count_lpdf;
  using stan::math::test::get_new_msg;
  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  std::vector<int> idata;
  std::vector<double> vlambda_d(1, lambda_d);

  double poisson_lpdf = stan::math::reduce_sum<count_lpdf<double>>(
      data, 5, get_new_msg(), vlambda_d, idata);
  double poisson_static_lpdf
      = stan::math::reduce_sum_static<count_lpdf<double>>(
          data, 5, get_new_msg(), vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);
  // NOTE:(Steve) This fails with EXPECT_DOUBLE_EQ at about 10e-7
  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);
  EXPECT_FLOAT_EQ(poisson_static_lpdf, poisson_lpdf_ref);
}

TEST(StanMathPrim_reduce_sum, grainsize) {
  using stan::math::test::get_new_msg;
  using stan::math::test::sum_lpdf;

  std::vector<int> data(5, 10);

  EXPECT_THROW(stan::math::reduce_sum<sum_lpdf>(data, 0, get_new_msg()),
               std::domain_error);

  EXPECT_THROW(stan::math::reduce_sum<sum_lpdf>(data, -1, get_new_msg()),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::reduce_sum<sum_lpdf>(data, 1, get_new_msg()));

  EXPECT_NO_THROW(
      stan::math::reduce_sum<sum_lpdf>(data, 3 * data.size(), get_new_msg()));
  EXPECT_THROW(stan::math::reduce_sum_static<sum_lpdf>(data, 0, get_new_msg()),
               std::domain_error);

  EXPECT_THROW(stan::math::reduce_sum_static<sum_lpdf>(data, -1, get_new_msg()),
               std::domain_error);

  EXPECT_NO_THROW(
      stan::math::reduce_sum_static<sum_lpdf>(data, 1, get_new_msg()));

  EXPECT_NO_THROW(stan::math::reduce_sum_static<sum_lpdf>(data, 3 * data.size(),
                                                          get_new_msg()));
}

TEST(StanMathPrim_reduce_sum, start_end_slice) {
  using stan::math::test::get_new_msg;
  using stan::math::test::start_end_lpdf;

  std::vector<int> data(5, 10);

  EXPECT_EQ(
      50, stan::math::reduce_sum<start_end_lpdf>(data, 1, get_new_msg(), data));
}

TEST(StanMathPrim_reduce_sum, std_vector_int_slice) {
  stan::math::test::test_slices(50, 10);
}

TEST(StanMathPrim_reduce_sum, std_vector_double_slice) {
  stan::math::test::test_slices(50.0, 10.0);
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_double_slice) {
  stan::math::test::test_slices(100.0, std::vector<double>(2, 10));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_int_slice) {
  stan::math::test::test_slices(100, std::vector<int>(2, 10));
}

TEST(StanMathPrim_reduce_sum, std_vector_eigen_vector_slice) {
  Eigen::VectorXd test_x = Eigen::VectorXd::Ones(2);
  stan::math::test::test_slices(10.0, test_x);
}

TEST(StanMathPrim_reduce_sum, std_vector_eigen_row_vector_slice) {
  Eigen::RowVectorXd test_x = Eigen::RowVectorXd::Ones(2);
  stan::math::test::test_slices(10.0, test_x);
}

TEST(StanMathPrim_reduce_sum, std_vector_eigen_matrix_slice) {
  Eigen::MatrixXd test_x = Eigen::MatrixXd::Ones(2, 2);
  stan::math::test::test_slices(20.0, test_x);
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_std_vector_int_slice) {
  stan::math::test::test_slices(
      200, std::vector<std::vector<int>>(2, std::vector<int>(2, 10.0)));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_std_vector_double_slice) {
  stan::math::test::test_slices(
      200.0, std::vector<std::vector<double>>(2, std::vector<double>(2, 10.0)));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_eigen_vector_slice) {
  stan::math::test::test_slices(
      20.0, std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(2)));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_eigen_row_vector_slice) {
  stan::math::test::test_slices(
      20.0, std::vector<Eigen::RowVectorXd>(2, Eigen::RowVectorXd::Ones(2)));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_eigen_matrix_slice) {
  stan::math::test::test_slices(
      40.0, std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Ones(2, 2)));
}

TEST(StanMathPrim_reduce_sum, no_args) {
  using stan::math::test::get_new_msg;
  using stan::math::test::sum_lpdf;

  std::vector<double> data(0);
  EXPECT_EQ(0.0, stan::math::reduce_sum_static<sum_lpdf>(
                     data, 1, stan::math::test::get_new_msg()))
      << "Failed for reduce_sum_static";
  EXPECT_EQ(0.0, stan::math::reduce_sum<sum_lpdf>(
                     data, 1, stan::math::test::get_new_msg()))
      << "Failed for reduce_sum";
}

TEST(StanMathPrim_reduce_sum, int_arg) {
  stan::math::test::test_slices(5 * (10 + 5), 10, 5);
}

TEST(StanMathPrim_reduce_sum, double_arg) {
  stan::math::test::test_slices(5 * (10.0 + 5.0), 10.0, 5.0);
}

TEST(StanMathPrim_reduce_sum, std_vector_int_arg) {
  stan::math::test::test_slices(5 * (10.0 + 5.0 * 10), 10.0,
                                std::vector<int>(5, 10));
}

TEST(StanMathPrim_reduce_sum, std_vector_double_arg) {
  stan::math::test::test_slices(5 * (10 + 5 * 10), 10.0,
                                std::vector<double>(5, 10));
}

TEST(StanMathPrim_reduce_sum, eigen_vector_arg) {
  stan::math::test::test_slices(5 * (10 + 5), 10.0, Eigen::VectorXd::Ones(5));
}

TEST(StanMathPrim_reduce_sum, eigen_row_vector_arg) {
  stan::math::test::test_slices(5 * (10 + 5), 10.0,
                                Eigen::RowVectorXd::Ones(5));
}

TEST(StanMathPrim_reduce_sum, eigen_matrix_arg) {
  stan::math::test::test_slices(5 * (10 + 5 * 5), 10.0,
                                Eigen::MatrixXd::Ones(5, 5));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_double_arg) {
  stan::math::test::test_slices(
      5 * (10 + 250), 10.0,
      std::vector<std::vector<double>>(5, std::vector<double>(5, 10.0)));
}

TEST(StanMathPrim_reduce_sum, std_vector_eigen_vector_arg) {
  stan::math::test::test_slices(
      5 * (10 + 10), 10.0,
      std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(5)));
}

TEST(StanMathPrim_reduce_sum, std_vector_eigen_row_vector_arg) {
  stan::math::test::test_slices(
      5 * (10 + 10), 10.0,
      std::vector<Eigen::RowVectorXd>(2, Eigen::RowVectorXd::Ones(5)));
}

TEST(StanMathPrim_reduce_sum, std_vector_eigen_matrix_arg) {
  stan::math::test::test_slices(
      5 * (10 + 50), 10.0,
      std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Ones(5, 5)));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_std_vector_double_arg) {
  stan::math::test::test_slices(5 * (10 + 1250), 10.0,
                                std::vector<std::vector<std::vector<double>>>(
                                    5, std::vector<std::vector<double>>(
                                           5, std::vector<double>(5, 10.0))));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_eigen_vector_arg) {
  stan::math::test::test_slices(
      5 * (10 + 20), 10.0,
      std::vector<std::vector<Eigen::VectorXd>>(
          2, std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(5))));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_eigen_row_vector_arg) {
  stan::math::test::test_slices(
      5 * (10 + 20), 10.0,
      std::vector<std::vector<Eigen::RowVectorXd>>(
          2, std::vector<Eigen::RowVectorXd>(2, Eigen::RowVectorXd::Ones(5))));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_eigen_matrix_arg) {
  stan::math::test::test_slices(
      5 * (10 + 60), 10.0,
      std::vector<std::vector<Eigen::MatrixXd>>(
          2, std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Ones(5, 3))));
}

TEST(StanMathPrim_reduce_sum, sum) {
  double answer = 5 + 5 * (1 + 1 + 5 + 5 + 5 + 5 + 25 + 10 + 10);
  stan::math::test::test_slices(
      answer, 1.0, 1, 1.0, std::vector<int>(5, 1), std::vector<double>(5, 1.0),
      Eigen::VectorXd::Ones(5), Eigen::RowVectorXd::Ones(5),
      Eigen::MatrixXd::Ones(5, 5),
      std::vector<std::vector<double>>(2, std::vector<double>(5, 1.0)),
      std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(5)));
}
