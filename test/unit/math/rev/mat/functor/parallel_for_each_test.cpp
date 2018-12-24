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
  using stan::math::ChainableStack;
  using stan::math::var;
  typedef boost::counting_iterator<int> count_iter;

  const int num_threads = stan::math::internal::get_num_threads();
  tbb::task_scheduler_init task_scheduler(num_threads);

  fun1 f;
  const int num_jobs = 1000;

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

  // next: copy the AD tape from the threads to the main AD tape
  typedef ChainableStack::AutodiffStackStorage local_ad_stack_t;

  local_ad_stack_t& ad_tape = ChainableStack::instance();
  std::for_each(ChainableStack::instance_.begin(),
                ChainableStack::instance_.end(),
                [&](local_ad_stack_t& local_instance) {
                  if (&local_instance == &ad_tape) {
                    std::cout << "Found main tape!" << std::endl;
                    return;
                  }
                  ad_tape.var_stack_.insert(ad_tape.var_stack_.end(),
                                            local_instance.var_stack_.begin(),
                                            local_instance.var_stack_.end());
                  local_instance.var_stack_.clear();
                  /*
                  it_t begin = local_instance.var_stack_.rbegin();
                  it_t end = local_instance.var_stack_.rend();
                  for (it_t it = begin; it != end; ++it) {
                    (*it)->chain();
                  }
                  */
                });

  for (int i = 0; i < num_jobs; ++i) {
    vector_d x_ref = x_ref_d[i];
    double fx_ref = value_of(fres[i]);
    vector<double> grad_fx_ref(2);
    // stan::math::set_zero_all_adjoints_global();
    // stan::math::grad_global(fres[i].vi_);
    stan::math::set_zero_all_adjoints();
    stan::math::grad(fres[i].vi_);
    grad_fx_ref[0] = x_ref_v[i](0).adj();
    grad_fx_ref[1] = x_ref_v[i](1).adj();

    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                    fx_ref);
    EXPECT_EQ(2, grad_fx_ref.size());
    EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref[0]);
    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref[1]);
  }
}
