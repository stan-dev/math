// the tests here check that map_rect_concurrent works correct as such
// we enforce that STAN_MPI is NOT defined

#ifdef STAN_MPI
#undef STAN_MPI
#endif

#include <gtest/gtest.h>
#include <stan/math/rev.hpp>

#include <test/unit/math/prim/functor/hard_work.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

#include <iostream>
#include <vector>

struct variadic_hard_work_1 {
  variadic_hard_work_1() {}
  template <typename T1>
  Eigen::Matrix<stan::return_type_t<T1>, Eigen::Dynamic, 1> operator()(
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
      std::ostream* msgs = 0) const {
    using result_type = stan::return_type_t<T1>;
    Eigen::Matrix<result_type, Eigen::Dynamic, 1> res;
    res.resize(2);
    res.setZero();
    res(0) = theta(0) * theta(0);
    res(1) = 2 * theta(1) * theta(0);
    return (res);
  }
};


STAN_REGISTER_MAP_RECT(0, hard_work)

struct map_rect : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  const std::size_t N = 100;

  virtual void SetUp() {
    set_n_threads(4);
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (std::size_t n = 0; n != N; ++n) {
      Eigen::VectorXd job_d(2);
      job_d << 0, n * n;
      job_params_d.push_back(job_d);
    }

    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
};

template <typename F_variadic, typename T_job_param, typename... Args>
struct rect_adapter {
  // mapping works as:
  // * all scalar vars are aggregated into the first shared argument
  // * all scalar datas are aggregated into the x_r and are thus
  // copied num_jobs times
  // * any non-scalar (vector, matrix, array) var is appended to the
  // shared argument in a flattened form
  // * any non-scalar double is appended to the x_r with copying
  // * any non-scalar int is appended to the x_i with copying

  using F_adapter_t = rect_adapter<F_variadic, T_job_param, Args...>;
  
  using return_t = Eigen::Matrix<stan::return_type_t<T_job_param, Args...>,
                                 Eigen::Dynamic, 1>;
  
  using job_param_t = std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>;
  using shared_param_t = Eigen::Matrix<stan::return_type_t<Args...>, Eigen::Dynamic, 1>;
  using job_data_double_t = std::vector<std::vector<double>>;
  using job_data_int_t = std::vector<std::vector<int>>;
  using rect_arg_t = std::tuple<shared_param_t, job_param_t, job_data_double_t, job_data_int_t>;

  // maps rectangular call to variadic F_variadic
  return_t operator()(const shared_param_t& shared_params,
                      const Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>&
                      job_param,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i,
                      std::ostream* msgs = nullptr) const {
    // take rectangular shape input and back-transform to variadic
    // form and call F_variadic with apply
    // TODO
    F_variadic f;
    return f(job_param);
  }

  // calls rectangular map_rect from variadic arguments and uses this
  // struct as adapter
  template<int call_id>
  return_t
  to_rect(const job_param_t& job_params, Args&&... args) const {
    const int num_jobs = job_params.size();
    shared_param_t shared;
    std::vector<double> x_r;
    std::vector<int> x_i;
    // do magic recursion which converts things "the right way"
    // TODO !!!
    return stan::math::map_rect<call_id, F_adapter_t>(shared, job_params, job_data_double_t(num_jobs, x_r), job_data_int_t(num_jobs, x_i));
  }
};

// it would be easier to only have to register
//rect_adapter<variadic_hard_work_1> as type.
//STAN_REGISTER_MAP_RECT(1, rect_adapter<variadic_hard_work_1, ...>)



template <int call_id, typename F_variadic,
          typename T_job_param, typename... Args>
Eigen::Matrix<stan::return_type_t<T_job_param, Args...>, Eigen::Dynamic, 1>
map_rect_variadic(
    const std::vector<Eigen::Matrix<T_job_param, Eigen::Dynamic, 1>>& job_params,
    std::ostream* msgs,
    Args&&... args) {

  using F_adapter_t = rect_adapter<F_variadic, T_job_param, Args...>;
  F_adapter_t adapter;
  
  return adapter.template to_rect<call_id>(job_params);
}



TEST_F(map_rect, variadic_1_eval_ok_d) {
  stan::math::vector_v shared_params_v = stan::math::to_var(shared_params_d);
  stan::math::vector_v res1 = map_rect_variadic<0, variadic_hard_work_1>(
      job_params_d, nullptr);

  for (std::size_t i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        stan::math::value_of(res1(j)),
        job_params_d[i](0) * job_params_d[i](0));
    EXPECT_FLOAT_EQ(stan::math::value_of(res1(j + 1)),
                    2 * job_params_d[i](1) * job_params_d[i](0)
                        );
  }
}


/*

TEST_F(map_rect, concurrent_eval_ok_dv) {
  std::vector<stan::math::vector_v> job_params_v;

  for (std::size_t i = 0; i < N; i++)
    job_params_v.push_back(stan::math::to_var(job_params_d[i]));

  stan::math::vector_v res1 = stan::math::map_rect<0, hard_work>(
      shared_params_d, job_params_v, x_r, x_i);

  for (std::size_t i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        stan::math::value_of(res1(j)),
        job_params_d[i](0) * job_params_d[i](0) + shared_params_d(0));
    EXPECT_FLOAT_EQ(stan::math::value_of(res1(j + 1)),
                    x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                        + 2 * shared_params_d(0) + shared_params_d(1));

    stan::math::set_zero_all_adjoints();
    res1(j).grad();
    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_, 2.0 * job_params_d[i](0));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_, 0.0);
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }

    stan::math::set_zero_all_adjoints();
    res1(j + 1).grad();
    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_,
                        x_r[i][0] * job_params_d[i](1));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_,
                        x_r[i][0] * job_params_d[i](0));
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }
  }
}

TEST_F(map_rect, concurrent_eval_ok_vv) {
  stan::math::vector_v shared_params_v = stan::math::to_var(shared_params_d);
  std::vector<stan::math::vector_v> job_params_v;

  for (std::size_t i = 0; i < N; i++)
    job_params_v.push_back(stan::math::to_var(job_params_d[i]));

  stan::math::vector_v res1 = stan::math::map_rect<0, hard_work>(
      shared_params_v, job_params_v, x_r, x_i);

  for (std::size_t i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        stan::math::value_of(res1(j)),
        job_params_d[i](0) * job_params_d[i](0) + shared_params_d(0));
    EXPECT_FLOAT_EQ(stan::math::value_of(res1(j + 1)),
                    x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                        + 2 * shared_params_d(0) + shared_params_d(1));

    stan::math::set_zero_all_adjoints();
    res1(j).grad();
    EXPECT_FLOAT_EQ(shared_params_v(0).vi_->adj_, 1.0);
    EXPECT_FLOAT_EQ(shared_params_v(1).vi_->adj_, 0.0);

    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_, 2.0 * job_params_d[i](0));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_, 0.0);
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }

    stan::math::set_zero_all_adjoints();
    res1(j + 1).grad();
    EXPECT_FLOAT_EQ(shared_params_v(0).vi_->adj_, 2.0);
    EXPECT_FLOAT_EQ(shared_params_v(1).vi_->adj_, 1.0);

    for (std::size_t k = 0; k < N; k++) {
      if (k == i) {
        EXPECT_FLOAT_EQ(job_params_v[i](0).vi_->adj_,
                        x_r[i][0] * job_params_d[i](1));
        EXPECT_FLOAT_EQ(job_params_v[i](1).vi_->adj_,
                        x_r[i][0] * job_params_d[i](0));
      } else {
        EXPECT_FLOAT_EQ(job_params_v[k](0).vi_->adj_, 0.0);
        EXPECT_FLOAT_EQ(job_params_v[k](1).vi_->adj_, 0.0);
      }
    }
  }
}
*/
