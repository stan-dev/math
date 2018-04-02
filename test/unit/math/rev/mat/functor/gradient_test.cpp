#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
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

TEST(AgradAutoDiff, gradient) {
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

// test threaded AD
TEST(AgradAutoDiff, gradient_threaded) {
  fun1 f;
  Matrix<double, Dynamic, 1> x(2);
  x << 5, 7;
  double fx_ref;
  Matrix<double, Dynamic, 1> grad_fx_ref;
  stan::math::gradient(f, x, fx_ref, grad_fx_ref);
  EXPECT_FLOAT_EQ(5 * 5 * 7 + 3 * 7 * 7, fx_ref);
  EXPECT_EQ(2, grad_fx_ref.size());
  EXPECT_FLOAT_EQ(2 * x(0) * x(1), grad_fx_ref(0));
  EXPECT_FLOAT_EQ(x(0) * x(0) + 3 * 2 * x(1), grad_fx_ref(1));

  auto thread_job = [&]() {
    double fx;
    VectorXd grad_fx;
    stan::math::gradient(fun1(), x, fx, grad_fx);
    VectorXd res(1+grad_fx.size());
    res(0) = fx;
    res.tail(grad_fx.size()) = grad_fx;
    return res;
  };

  std::vector<std::future<VectorXd>> ad_futures;

  for (std::size_t i=0; i < 100; i++) {
    // the use pattern in stan-math will be to defer the first job in
    // order to make the main thread to some work which is why we
    // alter the execution policy here
    ad_futures.emplace_back(std::async(i == 0 ? std::launch::deferred : std::launch::async,
                                       thread_job));
  }

  for (std::size_t i=0; i < 100; i++) {
    const VectorXd& ad_result = ad_futures[i].get();
    double fx_job = ad_result(0);
    VectorXd grad_fx_job = ad_result.tail(ad_result.size()-1);

    EXPECT_FLOAT_EQ(fx_ref, fx_job);
    EXPECT_EQ(grad_fx_ref.size(), grad_fx_job.size());
    EXPECT_FLOAT_EQ(grad_fx_ref(0), grad_fx_job(0));
    EXPECT_FLOAT_EQ(grad_fx_ref(1), grad_fx_job(1));
  }
}

stan::math::var sum_and_throw(const Matrix<stan::math::var, Dynamic, 1>& x) {
  stan::math::var y = 0;
  for (int i = 0; i < x.size(); ++i)
    y += x(i);
  throw std::domain_error("fooey");
  return y;
}

TEST(AgradAutoDiff, RecoverMemory) {
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
  EXPECT_LT(stan::math::ChainableStack::context().memalloc_.bytes_allocated(), 100000);
}
