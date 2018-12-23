#include <stan/math/rev/mat.hpp>
#include <stan/math/parallel/for_each.hpp>
#include <stan/math/parallel/get_num_threads.hpp>
#include <gtest/gtest.h>

#include <tbb/task_scheduler_init.h>

#include <boost/iterator/counting_iterator.hpp>

#include <stdexcept>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using stan::math::value_of;
using stan::math::vector_d;
using stan::math::vector_v;
using std::vector;

// fun1(x, y) = (x^2 * y) + (3 * y^2)
struct fun1 {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1);
  }
};

// test threaded AD if enabled
TEST(AgradAutoDiff, parallel_for_each) {
  using stan::math::var;
  typedef boost::counting_iterator<int> count_iter;

  const int num_threads = stan::math::internal::get_num_threads();
  tbb::task_scheduler_init task_scheduler(num_threads);

  fun1 f;
  const int num_jobs = 10000;

  vector<vector_v> x_ref_v(num_jobs);
  vector<vector_d> x_ref_d(num_jobs);

  for (int i = 0; i < num_jobs; ++i) {
    vector_v xi(2);
    xi << 5, 7 + i;
    x_ref_v[i] = xi;
    x_ref_d[i] = value_of(xi);
  }

  vector<var> fres(num_jobs);

  auto apply_f = [&](int i) -> void { fres[i] = f(x_ref_v[i]); };

#ifdef STAN_THREADS
  constexpr std_par::execution::parallel_unsequenced_policy exec_policy
      = std_par::execution::par_unseq;
#else
  constexpr std_par::execution::sequenced_policy exec_policy
      = std_par::execution::seq;
#endif

  std_par::for_each(exec_policy, count_iter(0), count_iter(num_jobs), apply_f);

  for (int i = 0; i < num_jobs; ++i) {
    vector_d x_ref = x_ref_d[i];
    double fx_ref = value_of(fres[i]);
    vector<double> grad_fx_ref(2);
    stan::math::set_zero_all_adjoints();
    fres[i].grad();
    grad_fx_ref[0] = x_ref_v[i](0).adj();
    grad_fx_ref[1] = x_ref_v[i](1).adj();

    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                    fx_ref);
    EXPECT_EQ(2, grad_fx_ref.size());
    EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref[0]);
    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref[1]);
  }

  /*
  auto thread_job = [&](double x1, double x2) {
    double fx;
    VectorXd x_local(2);
    x_local << x1, x2;
    VectorXd grad_fx;
    stan::math::gradient(fun1(), x_local, fx, grad_fx);
    VectorXd res(1 + grad_fx.size());
    res(0) = fx;
    res.tail(grad_fx.size()) = grad_fx;
    return res;
  };

  // schedule a bunch of jobs which all do the same
  std::vector<std::future<VectorXd>> ad_futures_ref;

  for (std::size_t i = 0; i < 100; i++) {
    // the use pattern in stan-math will be to defer the first job in
    // order to make the main thread do some work which is why we
    // alter the execution policy here
    ad_futures_ref.emplace_back(std::async(i == 0 ? std::launch::deferred
#ifndef STAN_THREADS
                                                  : std::launch::deferred,
#else
                                                  : std::launch::async,
#endif
                                           thread_job, x_ref(0), x_ref(1)));
  }

  // and schedule a bunch of jobs which all do different things (all
  // at the same time)
  std::vector<std::future<VectorXd>> ad_futures_local;

  for (std::size_t i = 0; i < 100; i++) {
    ad_futures_local.emplace_back(std::async(i == 0 ? std::launch::deferred
#ifndef STAN_THREADS
                                                    : std::launch::deferred,
#else
                                                    : std::launch::async,
#endif
                                             thread_job, 1.0 * i, 2.0 * i));
  }

  for (std::size_t i = 0; i < 100; i++) {
    const VectorXd& ad_result = ad_futures_ref[i].get();
    double fx_job = ad_result(0);
    VectorXd grad_fx_job = ad_result.tail(ad_result.size() - 1);

    EXPECT_FLOAT_EQ(fx_ref, fx_job);
    EXPECT_EQ(grad_fx_ref.size(), grad_fx_job.size());
    EXPECT_FLOAT_EQ(grad_fx_ref(0), grad_fx_job(0));
    EXPECT_FLOAT_EQ(grad_fx_ref(1), grad_fx_job(1));
  }

  for (std::size_t i = 0; i < 100; i++) {
    const VectorXd& ad_result = ad_futures_local[i].get();
    double fx_job = ad_result(0);
    VectorXd x_local(2);
    x_local << 1.0 * i, 2.0 * i;
    VectorXd grad_fx_job = ad_result.tail(ad_result.size() - 1);

    EXPECT_FLOAT_EQ(
        x_local(0) * x_local(0) * x_local(1) + 3 * x_local(1) * x_local(1),
        fx_job);
    EXPECT_EQ(2, grad_fx_job.size());
    EXPECT_FLOAT_EQ(2 * x_local(0) * x_local(1), grad_fx_job(0));
    EXPECT_FLOAT_EQ(x_local(0) * x_local(0) + 3 * 2 * x_local(1),
                    grad_fx_job(1));
  }
  */
}
