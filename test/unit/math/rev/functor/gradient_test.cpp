#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <stdexcept>
#include <vector>
#include <thread>
#include <future>

using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// fun1(x, y) = (x^2 * y) + (3 * y^2)
struct fun1 {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1);
  }
};

TEST(RevFunctor, gradient) {
  fun1 f;
  Matrix<double, Dynamic, 1> x(2);
  x << 5, 7;
  double fx;
  Matrix<double, Dynamic, 1> grad_fx;
  stan::math::gradient(f, x, fx, grad_fx);
  EXPECT_FLOAT_EQ(5 * 5 * 7 + 3 * 7 * 7, fx);
  EXPECT_EQ(2, grad_fx.size());
  EXPECT_FLOAT_EQ(2 * x(0) * x(1), grad_fx(0));
  EXPECT_FLOAT_EQ(x(0) * x(0) + 3 * 2 * x(1), grad_fx(1));
}

TEST(RevFunctor, gradient_input_rowvector) {
  fun1 f;
  Matrix<double, 1, Dynamic> x(2);
  x << 5, 7;
  double fx;
  std::vector<double> grad_fx(2, 0);
  stan::math::gradient(f, x, fx, std::begin(grad_fx), std::end(grad_fx));
  EXPECT_FLOAT_EQ(5 * 5 * 7 + 3 * 7 * 7, fx);
  EXPECT_EQ(2, grad_fx.size());
  EXPECT_FLOAT_EQ(2 * x(0) * x(1), grad_fx[0]);
  EXPECT_FLOAT_EQ(x(0) * x(0) + 3 * 2 * x(1), grad_fx[1]);
}

TEST(RevFunctor, gradient_array) {
  fun1 f;
  Matrix<double, Dynamic, 1> x(2);
  x << 5, 7;
  double fx;
  std::vector<double> grad_fx(2, 0);
  stan::math::gradient(f, x, fx, std::begin(grad_fx), std::end(grad_fx));
  EXPECT_FLOAT_EQ(5 * 5 * 7 + 3 * 7 * 7, fx);
  EXPECT_EQ(2, grad_fx.size());
  EXPECT_FLOAT_EQ(2 * x(0) * x(1), grad_fx[0]);
  EXPECT_FLOAT_EQ(x(0) * x(0) + 3 * 2 * x(1), grad_fx[1]);
}

// test threaded AD if enabled
TEST(RevFunctor, gradient_threaded) {
  fun1 f;
  Matrix<double, Dynamic, 1> x_ref(2);
  x_ref << 5, 7;
  double fx_ref;
  Matrix<double, Dynamic, 1> grad_fx_ref;
  stan::math::gradient(f, x_ref, fx_ref, grad_fx_ref);
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                  fx_ref);
  EXPECT_EQ(2, grad_fx_ref.size());
  EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref(0));
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref(1));

  auto thread_job = [&](double x1, double x2) {
    stan::math::ChainableStack thread_instance;
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
}

// test threaded AD if enabled
TEST(RevFunctor, gradient_array_threaded) {
  fun1 f;
  Matrix<double, Dynamic, 1> x_ref(2);
  x_ref << 5, 7;
  double fx_ref;
  std::vector<double> grad_fx_ref(2, 0);
  stan::math::gradient(f, x_ref, fx_ref, std::begin(grad_fx_ref),
                       std::end(grad_fx_ref));
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                  fx_ref);
  EXPECT_EQ(2, grad_fx_ref.size());
  EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref[0]);
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref[1]);

  auto thread_job = [&](double x1, double x2) {
    stan::math::ChainableStack thread_instance;
    double fx;
    VectorXd x_local(2);
    x_local << x1, x2;
    std::vector<double> grad_fx(2, 0);
    stan::math::gradient(fun1(), x_local, fx, std::begin(grad_fx),
                         std::end(grad_fx));
    VectorXd res(1 + x_local.size());
    res(0) = fx;
    for (size_t i = 0; i < x_local.size(); i++) {
      res(i + 1) = grad_fx[i];
    }
    return res;
  };

  // schedule a bunch of jobs which all do the same
  std::vector<std::future<VectorXd>> ad_futures_ref;

  for (std::size_t i = 0; i < 100; i++) {
    /*
     * the use pattern in stan-math will be to defer the first job in
     * order to make the main thread do some work which is why we
     * alter the execution policy here
     */
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
    EXPECT_FLOAT_EQ(grad_fx_ref[0], grad_fx_job(0));
    EXPECT_FLOAT_EQ(grad_fx_ref[1], grad_fx_job(1));
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
}

// test threaded AD through the Intel TBB whenever threading is used
#ifdef STAN_THREADS
TEST(RevFunctor, gradient_threaded_tbb) {
  fun1 f;
  Matrix<double, Dynamic, 1> x_ref(2);
  x_ref << 5, 7;
  double fx_ref;
  Matrix<double, Dynamic, 1> grad_fx_ref;
  stan::math::gradient(f, x_ref, fx_ref, grad_fx_ref);
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                  fx_ref);
  EXPECT_EQ(2, grad_fx_ref.size());
  EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref(0));
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref(1));

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
  std::vector<VectorXd> ad_ref(100);

  using tbb::blocked_range;

  tbb::parallel_for(blocked_range<std::size_t>(0, 100),
                    [&](const blocked_range<size_t>& r) {
                      for (std::size_t i = r.begin(); i != r.end(); ++i)
                        ad_ref[i] = thread_job(x_ref(0), x_ref(1));
                    });

  for (std::size_t i = 0; i < 100; i++) {
    const VectorXd& ad_result = ad_ref[i];
    double fx_job = ad_result(0);
    VectorXd grad_fx_job = ad_result.tail(ad_result.size() - 1);

    EXPECT_FLOAT_EQ(fx_ref, fx_job);
    EXPECT_EQ(grad_fx_ref.size(), grad_fx_job.size());
    EXPECT_FLOAT_EQ(grad_fx_ref(0), grad_fx_job(0));
    EXPECT_FLOAT_EQ(grad_fx_ref(1), grad_fx_job(1));
  }

  // and schedule a bunch of jobs which all do different things (all
  // at the same time)
  std::vector<VectorXd> ad_local(100);

  tbb::parallel_for(blocked_range<std::size_t>(0, 100),
                    [&](const blocked_range<size_t>& r) {
                      for (std::size_t i = r.begin(); i != r.end(); ++i)
                        ad_local[i] = thread_job(1.0 * i, 2.0 * i);
                    });

  for (std::size_t i = 0; i < 100; i++) {
    const VectorXd& ad_result = ad_local[i];
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
}

TEST(RevFunctor, gradient_array_threaded_tbb) {
  fun1 f;
  Matrix<double, Dynamic, 1> x_ref(2);
  x_ref << 5, 7;
  double fx_ref;
  std::vector<double> grad_fx_ref(2, 0);
  stan::math::gradient(f, x_ref, fx_ref, std::begin(grad_fx_ref),
                       std::end(grad_fx_ref));
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                  fx_ref);
  EXPECT_EQ(2, grad_fx_ref.size());
  EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref[0]);
  EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref[1]);

  auto thread_job = [&](double x1, double x2) {
    double fx;
    VectorXd x_local(2);
    x_local << x1, x2;
    VectorXd grad_fx;
    stan::math::gradient(fun1(), x_local, fx, grad_fx);
    VectorXd res(1 + x_local.size());
    res(0) = fx;
    for (size_t i = 0; i < x_local.size(); i++) {
      res(i + 1) = grad_fx[i];
    }
    return res;
  };

  // schedule a bunch of jobs which all do the same
  std::vector<VectorXd> ad_ref(100);

  using tbb::blocked_range;

  tbb::parallel_for(blocked_range<std::size_t>(0, 100),
                    [&](const blocked_range<size_t>& r) {
                      for (std::size_t i = r.begin(); i != r.end(); ++i)
                        ad_ref[i] = thread_job(x_ref(0), x_ref(1));
                    });

  for (std::size_t i = 0; i < 100; i++) {
    const VectorXd& ad_result = ad_ref[i];
    double fx_job = ad_result(0);
    VectorXd grad_fx_job = ad_result.tail(ad_result.size() - 1);

    EXPECT_FLOAT_EQ(fx_ref, fx_job);
    EXPECT_EQ(grad_fx_ref.size(), grad_fx_job.size());
    EXPECT_FLOAT_EQ(grad_fx_ref[0], grad_fx_job(0));
    EXPECT_FLOAT_EQ(grad_fx_ref[1], grad_fx_job(1));
  }

  // and schedule a bunch of jobs which all do different things (all
  // at the same time)
  std::vector<VectorXd> ad_local(100);

  tbb::parallel_for(blocked_range<std::size_t>(0, 100),
                    [&](const blocked_range<size_t>& r) {
                      for (std::size_t i = r.begin(); i != r.end(); ++i)
                        ad_local[i] = thread_job(1.0 * i, 2.0 * i);
                    });

  for (std::size_t i = 0; i < 100; i++) {
    const VectorXd& ad_result = ad_local[i];
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
}
#endif

stan::math::var sum_and_throw(const Matrix<stan::math::var, Dynamic, 1>& x) {
  stan::math::var y = 0;
  for (int i = 0; i < x.size(); ++i)
    y += x(i);
  throw std::domain_error("fooey");
}

TEST(RevFunctor, RecoverMemory) {
  using Eigen::VectorXd;
  for (int i = 0; i < 100000; ++i) {
    try {
      VectorXd x(5);
      x << 1, 2, 3, 4, 5;
      double fx;
      VectorXd grad_fx;
      stan::math::gradient(sum_and_throw, x, fx, grad_fx);
    } catch (const std::domain_error& e) {
      // ignore me
    }
  }
  // depends on starting allocation of 65K not being exceeded
  // without recovery_memory in autodiff::apply_recover(), takes 67M
  EXPECT_LT(stan::math::ChainableStack::instance_->memalloc_.bytes_allocated(),
            100000);
}

TEST(RevFunctor, RecoverMemory_gradient_array) {
  using Eigen::VectorXd;
  for (int i = 0; i < 100000; ++i) {
    try {
      VectorXd x(5);
      x << 1, 2, 3, 4, 5;
      double fx;
      std::vector<double> grad_fx(5, 0);
      stan::math::gradient(sum_and_throw, x, fx, std::begin(grad_fx),
                           std::end(grad_fx));
    } catch (const std::domain_error& e) {
      // ignore me
    }
  }
  // depends on starting allocation of 65K not being exceeded
  // without recovery_memory in autodiff::apply_recover(), takes 67M
  EXPECT_LT(stan::math::ChainableStack::instance_->memalloc_.bytes_allocated(),
            100000);
}

TEST(RevFunctor, gradientBoundaryConds) {
  VectorXd x(5);
  using stan::math::gradient;
  x << 1, 2, 3, 4, 5;
  double fx;
  std::vector<double> grad_fx(5, 0);
  EXPECT_NO_THROW(gradient([](const auto& x) { return stan::math::sum(x); }, x,
                           fx, std::begin(grad_fx), std::end(grad_fx)));
  EXPECT_THROW(gradient([](const auto& x) { return stan::math::sum(x); }, x, fx,
                        std::begin(grad_fx) + 1, std::end(grad_fx)),
               std::invalid_argument);
  EXPECT_THROW(gradient([](const auto& x) { return stan::math::sum(x); }, x, fx,
                        std::begin(grad_fx), std::end(grad_fx) + 1),
               std::invalid_argument);
}
