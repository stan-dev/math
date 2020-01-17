#include <stan/math/prim/core/init_threadpool_tbb.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

// reduce functor which is the BinaryFunction
// here we use iterators which represent integer indices
template <typename T>
struct count_lpdf {
  count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<int>& sub_slice,
                      const std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::poisson_lpmf(sub_slice, lambda[0]);
  }
};

TEST(StanMath_reduce_sum, value) {
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

  double poisson_lpdf
      = stan::math::reduce_sum<count_lpdf<double>>(data, 5, vlambda_d, idata);

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
                      const std::vector<int>& sub_slice,
                      const std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::reduce_sum<count_lpdf<T>>(sub_slice, 5, lambda, idata);
  }
};

TEST(StanMath_reduce_sum, nesting_value) {
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

  double poisson_lpdf
      = stan::math::reduce_sum<count_lpdf<double>>(data, 5, vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);
}

struct int_slice_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
                      const std::vector<int>& sub_slice) const {
    return stan::math::sum(sub_slice);
  }
};

struct double_slice_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
                      const std::vector<double>& sub_slice) const {
    return stan::math::sum(sub_slice);
  }
};

TEST(StanMath_reduce_sum, int_slice) {
  stan::math::init_threadpool_tbb();

  std::vector<int> data(5, 10);
  
  EXPECT_EQ(50, stan::math::reduce_sum<int_slice_lpdf>(data, 0));
}

TEST(StanMath_reduce_sum, double_slice) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  
  EXPECT_DOUBLE_EQ(50.0, stan::math::reduce_sum<double_slice_lpdf>(data, 0));
}

struct int_arg_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const int& arg) const {
    return stan::math::sum(sub_slice) + sub_slice.size() * arg;
  }
};

struct double_arg_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const double& arg) const {
    return stan::math::sum(sub_slice) + sub_slice.size() * arg;
  }
};

struct std_vector_int_arg_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const std::vector<int>& arg) const {
    return stan::math::sum(sub_slice) + sub_slice.size() * stan::math::sum(arg);
  }
};

struct std_vector_double_arg_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const std::vector<double>& arg) const {
    return stan::math::sum(sub_slice) + sub_slice.size() * stan::math::sum(arg);
  }
};

struct eigen_vector_arg_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const Eigen::VectorXd& arg) const {
    return stan::math::sum(sub_slice) + sub_slice.size() * stan::math::sum(arg);
  }
};

struct eigen_row_vector_arg_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const Eigen::RowVectorXd& arg) const {
    return stan::math::sum(sub_slice) + sub_slice.size() * stan::math::sum(arg);
  }
};

struct eigen_matrix_arg_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const Eigen::MatrixXd& arg) const {
    return stan::math::sum(sub_slice) + sub_slice.size() * stan::math::sum(arg);
  }
};

TEST(StanMath_reduce_sum, int_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  int arg = 5;
  
  EXPECT_DOUBLE_EQ(5 * (10 + 5), stan::math::reduce_sum<int_arg_lpdf>(data, 0, arg));
}

TEST(StanMath_reduce_sum, double_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  double arg = 5.0;
  
  EXPECT_DOUBLE_EQ(5 * (10.0 + 5.0), stan::math::reduce_sum<double_arg_lpdf>(data, 0, arg));
}

TEST(StanMath_reduce_sum, std_vector_int_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  std::vector<int> arg(5, 10);
  
  EXPECT_DOUBLE_EQ(5 * (10 + 5 * 10), stan::math::reduce_sum<std_vector_int_arg_lpdf>(data, 0, arg));
}

TEST(StanMath_reduce_sum, std_vector_double_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  std::vector<double> arg(5, 10.0);
  
  EXPECT_DOUBLE_EQ(5 * (10 + 5 * 10), stan::math::reduce_sum<std_vector_double_arg_lpdf>(data, 0, arg));
}

TEST(StanMath_reduce_sum, eigen_vector_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(5);
  
  EXPECT_DOUBLE_EQ(5 * (10 + 5), stan::math::reduce_sum<eigen_vector_arg_lpdf>(data, 0, arg));
}

TEST(StanMath_reduce_sum, eigen_row_vector_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(5);
  
  EXPECT_DOUBLE_EQ(5 * (10 + 5), stan::math::reduce_sum<eigen_row_vector_arg_lpdf>(data, 0, arg));
}

TEST(StanMath_reduce_sum, eigen_matrix_arg) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 10.0);
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(5, 5);
  
  EXPECT_DOUBLE_EQ(5 * (10 + 5 * 5), stan::math::reduce_sum<eigen_matrix_arg_lpdf>(data, 0, arg));
}

struct sum_lpdf {
  inline auto operator()(std::size_t start, std::size_t end,
			 const std::vector<double>& sub_slice,
			 const int& arg1,
			 const double& arg2,
			 const std::vector<int>& arg3,
			 const std::vector<double>& arg4,
			 const Eigen::VectorXd& arg5,
			 const Eigen::RowVectorXd& arg6,
			 const Eigen::MatrixXd& arg7) const {
    return stan::math::sum(sub_slice) +
      sub_slice.size() * (arg1 + arg2 +
			  stan::math::sum(arg3) + stan::math::sum(arg4) +
			  stan::math::sum(arg5) + stan::math::sum(arg6) +
			  stan::math::sum(arg7));
  }
};

TEST(StanMath_reduce_sum, sum) {
  stan::math::init_threadpool_tbb();

  std::vector<double> data(5, 1.0);
  int arg1 = 1;
  double arg2 = 1.0;
  std::vector<int> arg3(5, 1);
  std::vector<double> arg4(5, 1.0);
  Eigen::VectorXd arg5 = Eigen::VectorXd::Ones(5);
  Eigen::RowVectorXd arg6 = Eigen::RowVectorXd::Ones(5);
  Eigen::MatrixXd arg7 = Eigen::MatrixXd::Ones(5, 5);
  
  EXPECT_DOUBLE_EQ(5 + 5 * (1 + 1 + 5 + 5 + 5 + 5 + 25),
		   stan::math::reduce_sum<sum_lpdf>(data, 0, arg1, arg2,
						    arg3, arg4, arg5, arg6,
						    arg7));
}
