#include <stan/math/prim/core/init_threadpool_tbb.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

std::ostream* msgs = nullptr;

// reduce functor which is the BinaryFunction
// here we use iterators which represent integer indices
template <typename T>
struct count_lpdf {
  count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<int>& sub_slice, std::ostream* msgs,
                      const std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::poisson_lpmf(sub_slice, lambda[0]);
  }
};

TEST(StanMathPrim_reduce_sum, value) {
  stan::math::init_threadpool_tbb();

  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  std::vector<int> idata;
  std::vector<double> vlambda_d(1, lambda_d);

  typedef boost::counting_iterator<std::size_t> count_iter;

  double poisson_lpdf = stan::math::reduce_sum<count_lpdf<double>>(
      data, 5, msgs, vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);
}

// ********************************
// test if nested parallelism works
// ********************************

template <typename T>
struct nesting_count_lpdf {
  nesting_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<int>& sub_slice, std::ostream* msgs,
                      const std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::reduce_sum<count_lpdf<T>>(sub_slice, 5, msgs, lambda,
                                                 idata);
  }
};

TEST(StanMathPrim_reduce_sum, nesting_value) {
  stan::math::init_threadpool_tbb();

  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  std::vector<int> idata;
  std::vector<double> vlambda_d(1, lambda_d);

  typedef boost::counting_iterator<std::size_t> count_iter;

  double poisson_lpdf = stan::math::reduce_sum<count_lpdf<double>>(
      data, 5, msgs, vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);
}

template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
T sum_(T arg) {
  return arg;
}

template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
auto sum_(EigMat&& arg) {
  return stan::math::sum(arg);
}

template <typename Vec, stan::require_std_vector_t<Vec>* = nullptr>
auto sum_(Vec&& arg) {
  stan::scalar_type_t<Vec> sum = 0;
  for (size_t i = 0; i < arg.size(); ++i) {
    sum += sum_(arg[i]);
  }
  return sum;
}

struct sum_lpdf {
  template <typename T, typename... Args>
  inline auto operator()(std::size_t start, std::size_t end, T&& sub_slice,
                         std::ostream* msgs, Args&&... args) const {
    using return_type = stan::return_type_t<T, Args...>;

    return sum_(sub_slice)
           + sub_slice.size()
                 * stan::math::sum(std::vector<return_type>{
                       return_type(sum_(std::forward<Args>(args)))...});
  }
};

TEST(StanMathPrim_reduce_sum, int_slice) {
  stan::math::init_threadpool_tbb();

  std::vector<int> data(5, 10);

  EXPECT_EQ(50, stan::math::reduce_sum<sum_lpdf>(data, 0, msgs));
}

TEST(StanMathPrim_reduce_sum, double_slice) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);

  EXPECT_DOUBLE_EQ(50.0,
                   stan::math::reduce_sum<sum_lpdf>(data, 0, msgs));
}

TEST(StanMathPrim_reduce_sum, std_vector_double_slice) {
  stan::math::init_threadpool_tbb();

  std::vector<std::vector<double>> data(5, std::vector<double>(2, 10.0));

  EXPECT_DOUBLE_EQ(100.0,
                   stan::math::reduce_sum<sum_lpdf>(data, 0, msgs));
}

struct start_end_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
                         const std::vector<int>&, std::ostream* msgs,
                         const std::vector<int>& data) const {
    int sum = 0;
    EXPECT_GE(start, 0);
    EXPECT_LE(end, data.size() - 1);
    for (size_t i = start; i <= end; i++) {
      sum += data[i];
    }
    return sum;
  }
};

TEST(StanMathPrim_reduce_sum, start_end_slice) {
  stan::math::init_threadpool_tbb();

  std::vector<int> data(5, 10);

  EXPECT_EQ(50, stan::math::reduce_sum<start_end_lpdf>(data, 0, msgs, data));
}

TEST(StanMathPrim_reduce_sum, int_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  int arg = 5;

  EXPECT_DOUBLE_EQ(5 * (10 + 5),
                   stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, double_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  double arg = 5.0;

  EXPECT_DOUBLE_EQ(5 * (10.0 + 5.0),
                   stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, std_vector_int_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  std::vector<int> arg(5, 10);

  EXPECT_DOUBLE_EQ(
      5 * (10 + 5 * 10),
      stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, std_vector_double_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  std::vector<double> arg(5, 10.0);

  EXPECT_DOUBLE_EQ(
      5 * (10 + 5 * 10),
      stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, eigen_vector_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(5);

  EXPECT_DOUBLE_EQ(5 * (10 + 5), stan::math::reduce_sum<sum_lpdf>(
                                     data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, eigen_row_vector_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(5);

  EXPECT_DOUBLE_EQ(
      5 * (10 + 5),
      stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, eigen_matrix_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(5, 5);

  EXPECT_DOUBLE_EQ(
      5 * (10 + 5 * 5),
      stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, std_vector_std_vector_double_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  std::vector<std::vector<double>> arg(5, std::vector<double>(5, 10.0));

  EXPECT_DOUBLE_EQ(
      5 * (10 + 250),
      stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, std_vector_eigen_vector_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(5));

  EXPECT_DOUBLE_EQ(5 * (10 + 10),
                   stan::math::reduce_sum<sum_lpdf>(
                       data, 0, msgs, arg));
}

TEST(StanMathPrim_reduce_sum, sum) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 1.0);
  int arg1 = 1;
  double arg2 = 1.0;
  std::vector<int> arg3(5, 1);
  std::vector<double> arg4(5, 1.0);
  Eigen::VectorXd arg5 = Eigen::VectorXd::Ones(5);
  Eigen::RowVectorXd arg6 = Eigen::RowVectorXd::Ones(5);
  Eigen::MatrixXd arg7 = Eigen::MatrixXd::Ones(5, 5);
  std::vector<std::vector<double>> arg8(2, arg4);
  std::vector<Eigen::VectorXd> arg9(2, arg5);

  EXPECT_DOUBLE_EQ(
      5 + 5 * (1 + 1 + 5 + 5 + 5 + 5 + 25 + 10 + 10),
      stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, arg1, arg2, arg3, arg4,
					arg5, arg6, arg7, arg8, arg9));
}
