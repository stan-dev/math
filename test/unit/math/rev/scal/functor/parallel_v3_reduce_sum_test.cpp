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
                      std::vector<int>& sub_slice, std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::poisson_lpmf(sub_slice, lambda[0]);
  }
};

TEST(parallel_v3_reduce_sum, value) {
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

  double poisson_lpdf = stan::math::parallel_reduce_sum<count_lpdf<double>>(
      data.begin(), data.end(), 0.0, 5, vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);

  std::cout << "ref value of poisson lpdf : " << poisson_lpdf_ref << std::endl;
  ;

  std::cout << "value of poisson lpdf : " << poisson_lpdf << std::endl;
  ;
}

TEST(parallel_v3_reduce_sum, gradient) {
  stan::math::init_threadpool_tbb();

  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  var lambda_v = lambda_d;

  std::vector<int> idata;
  std::vector<var> vlambda_v(1, lambda_v);

  var poisson_lpdf = stan::math::parallel_reduce_sum<count_lpdf<var>>(
      data.begin(), data.end(), var(0.0), 5, vlambda_v, idata);

  var lambda_ref = lambda_d;
  var poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_ref);

  EXPECT_FLOAT_EQ(value_of(poisson_lpdf), value_of(poisson_lpdf_ref));

  stan::math::grad(poisson_lpdf_ref.vi_);
  const double lambda_ref_adj = lambda_ref.adj();

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf.vi_);
  const double lambda_adj = lambda_v.adj();

  EXPECT_FLOAT_EQ(lambda_adj, lambda_ref_adj);

  std::cout << "ref value of poisson lpdf : " << poisson_lpdf_ref.val()
            << std::endl;
  ;
  std::cout << "ref gradient wrt to lambda: " << lambda_ref_adj << std::endl;
  ;

  std::cout << "value of poisson lpdf : " << poisson_lpdf.val() << std::endl;
  ;
  std::cout << "gradient wrt to lambda: " << lambda_adj << std::endl;
  ;
  stan::math::recover_memory();
}

// ********************************
// test if nested parallelism works
// ********************************

template <typename T>
struct nesting_count_lpdf {
  nesting_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      std::vector<int>& sub_slice, std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::parallel_reduce_sum<count_lpdf<T>>(
        sub_slice.begin(), sub_slice.end(), T(0.0), 5, lambda, idata);
  }
};

TEST(parallel_v3_reduce_sum, nesting_value) {
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

  double poisson_lpdf = stan::math::parallel_reduce_sum<count_lpdf<double>>(
      data.begin(), data.end(), 0.0, 5, vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);

  std::cout << "ref value of poisson lpdf : " << poisson_lpdf_ref << std::endl;
  ;

  std::cout << "value of poisson lpdf : " << poisson_lpdf << std::endl;
  ;
}

TEST(parallel_v3_reduce_sum, nesting_gradient) {
  stan::math::init_threadpool_tbb();

  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  var lambda_v = lambda_d;

  std::vector<int> idata;
  std::vector<var> vlambda_v(1, lambda_v);

  var poisson_lpdf = stan::math::parallel_reduce_sum<nesting_count_lpdf<var>>(
      data.begin(), data.end(), var(0.0), 5, vlambda_v, idata);

  var lambda_ref = lambda_d;
  var poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_ref);

  EXPECT_FLOAT_EQ(value_of(poisson_lpdf), value_of(poisson_lpdf_ref));

  stan::math::grad(poisson_lpdf_ref.vi_);
  const double lambda_ref_adj = lambda_ref.adj();

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf.vi_);
  const double lambda_adj = lambda_v.adj();

  EXPECT_FLOAT_EQ(lambda_adj, lambda_ref_adj);

  std::cout << "ref value of poisson lpdf : " << poisson_lpdf_ref.val()
            << std::endl;
  ;
  std::cout << "ref gradient wrt to lambda: " << lambda_ref_adj << std::endl;
  ;

  std::cout << "value of poisson lpdf : " << poisson_lpdf.val() << std::endl;
  ;
  std::cout << "gradient wrt to lambda: " << lambda_adj << std::endl;
  ;

  stan::math::recover_memory();
}
