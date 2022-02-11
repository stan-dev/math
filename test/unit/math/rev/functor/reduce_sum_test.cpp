#include <stan/math.hpp>
#include <test/unit/math/prim/functor/reduce_sum_util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>
#include <set>

TEST(StanMathRev_reduce_sum, no_args) {
  using stan::math::var;
  using stan::math::test::get_new_msg;
  using stan::math::test::sum_lpdf;
  std::vector<var> data(0);
  EXPECT_EQ(0.0, stan::math::reduce_sum_static<sum_lpdf>(
                     data, 1, stan::math::test::get_new_msg())
                     .val())
      << "Failed for reduce_sum_static";
  EXPECT_EQ(0.0, stan::math::reduce_sum<sum_lpdf>(
                     data, 1, stan::math::test::get_new_msg())
                     .val())
      << "Failed for reduce_sum";
}

TEST(StanMathRev_reduce_sum, value) {
  using stan::math::test::count_lpdf;
  using stan::math::test::get_new_msg;
  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  std::vector<int> idata;
  std::vector<double> vlambda_d(1, lambda_d);

  double poisson_lpdf = stan::math::reduce_sum<count_lpdf<double>>(
      data, 5, get_new_msg(), vlambda_d, idata);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_d);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref)
      << "ref value of poisson lpdf : " << poisson_lpdf_ref << std::endl
      << "value of poisson lpdf : " << poisson_lpdf << std::endl;

  double poisson_lpdf_static
      = stan::math::reduce_sum_static<count_lpdf<double>>(
          data, 5, get_new_msg(), vlambda_d, idata);

  EXPECT_FLOAT_EQ(poisson_lpdf_static, poisson_lpdf_ref);
}

TEST(StanMathRev_reduce_sum, gradient) {
  using stan::math::var;
  using stan::math::test::count_lpdf;
  using stan::math::test::get_new_msg;

  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  var lambda_v = lambda_d;

  std::vector<int> idata;
  std::vector<var> vlambda_v(1, lambda_v);

  var poisson_lpdf = stan::math::reduce_sum<count_lpdf<var>>(
      data, 5, get_new_msg(), vlambda_v, idata);

  var lambda_ref = lambda_d;
  var poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_ref);

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

  var poisson_lpdf_static = stan::math::reduce_sum_static<count_lpdf<var>>(
      data, 5, get_new_msg(), vlambda_v, idata);

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf_static.vi_);
  const double lambda_adj_static = lambda_v.adj();
  EXPECT_FLOAT_EQ(lambda_adj_static, lambda_ref_adj);
  stan::math::recover_memory();
}

TEST(StanMathRev_reduce_sum, grainsize) {
  using stan::math::var;
  using stan::math::test::count_lpdf;
  using stan::math::test::get_new_msg;

  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  var lambda_v = lambda_d;

  std::vector<int> idata;
  std::vector<var> vlambda_v(1, lambda_v);

  EXPECT_THROW(stan::math::reduce_sum<count_lpdf<var>>(data, 0, get_new_msg(),
                                                       vlambda_v, idata),
               std::domain_error);

  EXPECT_THROW(stan::math::reduce_sum<count_lpdf<var>>(data, -1, get_new_msg(),
                                                       vlambda_v, idata),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::reduce_sum<count_lpdf<var>>(
      data, 1, get_new_msg(), vlambda_v, idata));

  EXPECT_NO_THROW(stan::math::reduce_sum<count_lpdf<var>>(
      data, 2 * elems, get_new_msg(), vlambda_v, idata));

  EXPECT_THROW(stan::math::reduce_sum_static<count_lpdf<var>>(
                   data, 0, get_new_msg(), vlambda_v, idata),
               std::domain_error);

  EXPECT_THROW(stan::math::reduce_sum_static<count_lpdf<var>>(
                   data, -1, get_new_msg(), vlambda_v, idata),
               std::domain_error);

  EXPECT_NO_THROW(stan::math::reduce_sum_static<count_lpdf<var>>(
      data, 1, get_new_msg(), vlambda_v, idata));

  EXPECT_NO_THROW(stan::math::reduce_sum_static<count_lpdf<var>>(
      data, 2 * elems, get_new_msg(), vlambda_v, idata));

  stan::math::recover_memory();
}

TEST(StanMathRev_reduce_sum, nesting_gradient) {
  using stan::math::var;
  using stan::math::test::get_new_msg;
  using stan::math::test::nesting_count_lpdf;

  double lambda_d = 10.0;
  const std::size_t elems = 10000;
  std::vector<int> data(elems);

  for (std::size_t i = 0; i != elems; ++i)
    data[i] = i;

  var lambda_v = lambda_d;

  std::vector<int> idata;
  std::vector<var> vlambda_v(1, lambda_v);

  var poisson_lpdf = stan::math::reduce_sum<nesting_count_lpdf<var>>(
      data, 5, get_new_msg(), vlambda_v, idata);

  var lambda_ref = lambda_d;
  var poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_ref);

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
      = stan::math::reduce_sum_static<nesting_count_lpdf<var>>(
          data, 5, get_new_msg(), vlambda_v, idata);

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf_static.vi_);
  const double lambda_adj_static = lambda_v.adj();

  EXPECT_FLOAT_EQ(lambda_adj_static, lambda_ref_adj);
  stan::math::recover_memory();
}

TEST(StanMathRev_reduce_sum, grouped_gradient) {
  using stan::math::var;
  using stan::math::test::get_new_msg;
  using stan::math::test::grouped_count_lpdf;

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

  var poisson_lpdf = stan::math::reduce_sum<grouped_count_lpdf<var>>(
      data, 5, get_new_msg(), vlambda_v, gidx);

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
      = stan::math::reduce_sum_static<grouped_count_lpdf<var>>(
          data, 5, get_new_msg(), vlambda_v, gidx);

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf_static.vi_);
  const double lambda_adj_static = lambda_v.adj();
  EXPECT_FLOAT_EQ(lambda_adj_static, lambda_ref_adj);
  stan::math::recover_memory();
}

TEST(StanMathRev_reduce_sum, grouped_gradient_eigen) {
  using stan::math::var;
  using stan::math::test::get_new_msg;
  using stan::math::test::grouped_count_lpdf;

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

  Eigen::Matrix<var, -1, 1> vlambda_v(groups);

  for (std::size_t i = 0; i != groups; ++i)
    vlambda_v[i] = i + 0.2;
  var lambda_v = vlambda_v[0];

  var poisson_lpdf = stan::math::reduce_sum<grouped_count_lpdf<var>>(
      data, 5, get_new_msg(), vlambda_v, gidx);

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
      = stan::math::reduce_sum_static<grouped_count_lpdf<var>>(
          data, 5, get_new_msg(), vlambda_v, gidx);

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf_static.vi_);
  const double lambda_adj_static = lambda_v.adj();

  EXPECT_FLOAT_EQ(lambda_adj_static, lambda_ref_adj);

  stan::math::recover_memory();
}

TEST(StanMathRev_reduce_sum, slice_group_gradient) {
  using stan::math::var;
  using stan::math::test::get_new_msg;
  using stan::math::test::slice_group_count_lpdf;

  double lambda_d = 10.0;
  const std::size_t groups = 10;
  const std::size_t elems_per_group = 1000;
  const std::size_t elems = groups * elems_per_group;

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

  std::vector<var> vlambda_v;

  for (std::size_t i = 0; i != groups; ++i)
    vlambda_v.push_back(i + 0.2);

  var lambda_v = vlambda_v[0];

  var poisson_lpdf = stan::math::reduce_sum<slice_group_count_lpdf<var>>(
      vlambda_v, 5, get_new_msg(), data, gsidx);

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
      = stan::math::reduce_sum_static<slice_group_count_lpdf<var>>(
          vlambda_v, 5, get_new_msg(), data, gsidx);

  stan::math::set_zero_all_adjoints();
  stan::math::grad(poisson_lpdf_static.vi_);
  const double lambda_adj_static = lambda_v.adj();

  EXPECT_FLOAT_EQ(lambda_adj_static, lambda_ref_adj);

  stan::math::recover_memory();
}

// reduce_sum needs still work fine whenever multiple arguments point
// to the same vari internally
TEST(StanMathRev_reduce_sum, linked_args) {
  using stan::math::test::arg_start_end_lpdf;
  using stan::math::test::get_new_msg;
  auto start_end_rs = [](auto&& arg1, auto&& arg2) {
    return stan::math::reduce_sum_static<arg_start_end_lpdf>(
        arg1, 1, get_new_msg(), arg2, arg2);
  };

  auto start_end = [](auto&& arg1, auto&& arg2) {
    arg_start_end_lpdf worker;
    return worker(arg1, 0, arg1.size() - 1, get_new_msg(), arg2, arg2);
  };

  std::vector<double> data(5, 1.0);

  stan::math::set_zero_all_adjoints();

  std::vector<stan::math::var> param2(data.begin(), data.end());
  std::vector<std::vector<stan::math::var>> param1(5, param2);

  stan::math::var res = start_end(param1, param2);

  stan::math::grad(res.vi_);

  double param2_adj = param2[0].adj();

  stan::math::set_zero_all_adjoints();

  std::vector<stan::math::var> param2_rs(data.begin(), data.end());
  std::vector<std::vector<stan::math::var>> param1_rs(5, param2_rs);

  stan::math::var res_rs = start_end_rs(param1_rs, param2_rs);

  stan::math::grad(res_rs.vi_);

  EXPECT_FLOAT_EQ(param2_rs[0].adj(), param2_adj);

  stan::math::recover_memory();
}
