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
                      std::vector<int> sub_slice, T& lambda) const {
    return stan::math::poisson_lpmf(sub_slice, lambda);
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

  typedef boost::counting_iterator<std::size_t> count_iter;

  double poisson_lpdf = stan::math::parallel_reduce_sum<count_lpdf<double>>(
      data.begin(), data.end(), 0.0, 5, lambda_d);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);
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

  var poisson_lpdf = stan::math::parallel_reduce_sum<count_lpdf<var>>(
      data.begin(), data.end(), var(0.0), 5, lambda_v);

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
}
