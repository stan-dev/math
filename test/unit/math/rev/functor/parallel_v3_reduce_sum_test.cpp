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

TEST(v3_reduce_sum, value) {
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

  double poisson_lpdf = stan::math::reduce_sum<count_lpdf<double>>(data, 5, vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);

  std::cout << "ref value of poisson lpdf : " << poisson_lpdf_ref << std::endl;

  std::cout << "value of poisson lpdf : " << poisson_lpdf << std::endl;
}

TEST(v3_reduce_sum, gradient) {
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

  var poisson_lpdf = stan::math::reduce_sum<count_lpdf<var>>(data, 5,
								 vlambda_v, idata);

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

  std::cout << "ref gradient wrt to lambda: " << lambda_ref_adj << std::endl;

  std::cout << "value of poisson lpdf : " << poisson_lpdf.val() << std::endl;

  std::cout << "gradient wrt to lambda: " << lambda_adj << std::endl;

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
                      const std::vector<int>& sub_slice,
                      const std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::reduce_sum<count_lpdf<T>>(sub_slice, 5, lambda, idata);
  }
};

TEST(v3_reduce_sum, nesting_gradient) {
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

  var poisson_lpdf = stan::math::reduce_sum<nesting_count_lpdf<var>>(
      data, 5, vlambda_v, idata);

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

  std::cout << "ref gradient wrt to lambda: " << lambda_ref_adj << std::endl;

  std::cout << "value of poisson lpdf : " << poisson_lpdf.val() << std::endl;

  std::cout << "gradient wrt to lambda: " << lambda_adj << std::endl;

  stan::math::recover_memory();
}

// ********************************
// basic performance test for a hierarchical model
// ********************************

template <typename T>
struct grouped_count_lpdf {
  grouped_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<int>& sub_slice,
                      const std::vector<T>& lambda,
                      const std::vector<int>& gidx) const {
    const std::size_t num_terms = end - start + 1;
    // std::cout << "sub-slice " << start << " - " << end << "; num_terms = " <<
    // num_terms << "; size = " << sub_slice.size() << std::endl;
    std::vector<T> lambda_slice(num_terms, 0.0);
    for (std::size_t i = 0; i != num_terms; ++i)
      lambda_slice[i] = lambda[gidx[start + i]];
    return stan::math::poisson_lpmf(sub_slice, lambda_slice);
  }
};

TEST(v3_reduce_sum, grouped_gradient) {
  stan::math::init_threadpool_tbb();

  double lambda_d = 10.0;
  const std::size_t groups = 10;
  const std::size_t elems_per_group = 1000;
  const std::size_t elems = groups * elems_per_group;
  const std::size_t num_iter = 1000;

  std::vector<int> data(elems);
  std::vector<int> gidx(elems);

  for (std::size_t i = 0; i != elems; ++i) {
    data[i] = i;
    gidx[i] = i / elems_per_group;
  }

  using stan::math::var;

  std::vector<var> vlambda_v;

  for (std::size_t i = 0; i != groups; ++i)
    vlambda_v.push_back(i + 0.2);

  var lambda_v = vlambda_v[0];

  var poisson_lpdf = stan::math::reduce_sum<grouped_count_lpdf<var>>(
      data, 5, vlambda_v, gidx);

  std::vector<var> vref_lambda_v;
  for (std::size_t i = 0; i != elems; ++i) {
    vref_lambda_v.push_back(vlambda_v[gidx[i]]);
  }
  var lambda_ref = vlambda_v[0];

  var poisson_lpdf_ref = stan::math::poisson_lpmf(data, vref_lambda_v);

  EXPECT_FLOAT_EQ(value_of(poisson_lpdf), value_of(poisson_lpdf_ref));

  stan::math::grad(poisson_lpdf_ref.vi_);
  const double lambda_ref_adj = lambda_ref.adj();

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf.vi_);
  const double lambda_adj = lambda_v.adj();

  EXPECT_FLOAT_EQ(lambda_adj, lambda_ref_adj);

  std::cout << "ref value of poisson lpdf : " << poisson_lpdf_ref.val()
            << std::endl;

  std::cout << "ref gradient wrt to lambda: " << lambda_ref_adj << std::endl;

  std::cout << "value of poisson lpdf : " << poisson_lpdf.val() << std::endl;

  std::cout << "gradient wrt to lambda: " << lambda_adj << std::endl;

  stan::math::recover_memory();
}

/*
// ********************************
// slice over the grouping variable which is a var
// ********************************

template <typename T>
struct slice_group_count_lpdf {
  slice_group_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<T>& lambda_slice,
                      const std::vector<int>& y,
                      const std::vector<int>& gsidx) const {
    const std::size_t num_groups = end - start + 1;
    T result = 0.0;
    for (std::size_t i = 0; i != num_groups; ++i) {
      std::vector<int> y_group(y.begin() + gsidx[start + i],
                               y.begin() + gsidx[start + i + 1]);
      result += stan::math::poisson_lpmf(y_group, lambda_slice[i]);
    }
    return result;
  }
};

TEST(v3_reduce_sum, slice_group_gradient) {
  stan::math::init_threadpool_tbb();

  double lambda_d = 10.0;
  const std::size_t groups = 10;
  const std::size_t elems_per_group = 1000;
  const std::size_t elems = groups * elems_per_group;
  const std::size_t num_iter = 1000;

  std::vector<int> data(elems);
  std::vector<int> gidx(elems);
  std::vector<int> gsidx(groups + 1);


  for (std::size_t i = 0, k = 0; i != groups; ++i) {
    gsidx[i] = k;
    for (std::size_t j = 0; j != elems_per_group; ++j, ++k) {
      data[k] = k;
      gidx[k] = i;
    }
    gsidx[i + 1] = k;
  }

  using stan::math::var;

  std::vector<var> vlambda_v;

  for (std::size_t i = 0; i != groups; ++i)
    vlambda_v.push_back(i + 0.2);

  var lambda_v = vlambda_v[0];

  var poisson_lpdf = stan::math::reduce_sum<slice_group_count_lpdf<var>>(
      vlambda_v, 5, data, gsidx);

  std::vector<var> vref_lambda_v;
  for (std::size_t i = 0; i != elems; ++i) {
    vref_lambda_v.push_back(vlambda_v[gidx[i]]);
  }
  var lambda_ref = vlambda_v[0];

  var poisson_lpdf_ref = stan::math::poisson_lpmf(data, vref_lambda_v);

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
*/
