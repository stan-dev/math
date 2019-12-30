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
                      std::vector<int> sub_slice) const {
    T lambda = 10.0;
    return stan::math::poisson_lpmf(sub_slice, lambda);
  }
};

TEST(parallel_v3_reduce_sum, gradient) {
  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  typedef boost::counting_iterator<std::size_t> count_iter;

  double poisson_lpdf = stan::math::parallel_reduce_sum<count_lpdf<double>>(
      data.begin(), data.end(), 0.0);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);
}

/*
TEST(parallel_v3_reduce_sum, gradient) {
  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  count_lpdf<double> reduce_op;

  var poisson_lpdf = stan::math::parallel_reduce_sum(
      count_iter(0), count_iter(elems), var(0.0), reduce_op);

  var lambda_ref = lambda_d;
  var poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_ref);

  EXPECT_FLOAT_EQ(value_of(poisson_lpdf), value_of(poisson_lpdf_ref));

  stan::math::grad(poisson_lpdf.vi_);
  stan::math::grad(poisson_lpdf_ref.vi_);

  EXPECT_FLOAT_EQ(poisson_lpdf.adj(), poisson_lpdf_ref.adj());
}
*/
