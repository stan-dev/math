#include <stan/math.hpp>
#include <test/unit/math/prim/functor/reduce_sum_util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>
#include <set>

struct closure_adapter {
  template<typename F, typename T_slice, typename... Args>
  auto operator()(const T_slice& subslice, std::size_t start,
                  std::size_t end, std::ostream* msgs,
                  const F& f, Args... args) {
    return f(msgs, subslice, start, end, args...);
  }
};

TEST(StanMathRev_reduce_sum, grouped_gradient_closure) {
  using stan::math::var;
  using stan::math::from_lambda;
  using stan::math::test::get_new_msg;

  double lambda_d = 10.0;
  const std::size_t groups = 10;
  const std::size_t elems_per_group = 1000;
  const std::size_t elems = groups * elems_per_group;

  std::vector<int> data(elems);
  std::vector<int> gidx(elems);

  for (std::size_t i = 0; i != elems; ++i) {
    data[i] = i;
    gidx[i] = i / elems_per_group;
  }

  std::vector<var> vlambda_v;

  for (std::size_t i = 0; i != groups; ++i)
    vlambda_v.push_back(i + 0.2);

  var lambda_v = vlambda_v[0];

  auto functor = from_lambda(
    [](auto& lambda, auto& slice, std::size_t start, std::size_t end, auto& gidx, std::ostream * msgs) {
      const std::size_t num_terms = end - start + 1;
      std::decay_t<decltype(lambda)> lambda_slice(num_terms);
      for (std::size_t i = 0; i != num_terms; ++i)
        lambda_slice[i] = lambda[gidx[start + i]];
      return stan::math::poisson_lpmf(slice, lambda_slice);
    }, vlambda_v);

  var poisson_lpdf = stan::math::reduce_sum(
      data, 5, get_new_msg(), functor, gidx);

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

  EXPECT_FLOAT_EQ(lambda_adj, lambda_ref_adj)
      << "ref value of poisson lpdf : " << poisson_lpdf_ref.val() << std::endl
      << "ref gradient wrt to lambda: " << lambda_ref_adj << std::endl
      << "value of poisson lpdf : " << poisson_lpdf.val() << std::endl
      << "gradient wrt to lambda: " << lambda_adj << std::endl;

  var poisson_lpdf_static
      = stan::math::reduce_sum_static(data, 5, get_new_msg(), functor, gidx);

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf_static.vi_);
  const double lambda_adj_static = lambda_v.adj();
  EXPECT_FLOAT_EQ(lambda_adj_static, lambda_ref_adj);
  stan::math::recover_memory();

  stan::math::recover_memory();
}
