#include <stan/math/rev/mat.hpp>
#include <stan/math/parallel/for_each.hpp>
#include <stan/math/parallel/get_num_threads.hpp>
#include <stan/math/rev/core/nest_chainablestack.hpp>
#include <stan/math/prim/scal/functor/parallel_reduce_sum.hpp>
#include <stan/math/rev/scal/functor/parallel_reduce_sum.hpp>
#include <gtest/gtest.h>

#define TBB_PREVIEW_LOCAL_OBSERVER 1

// #include <tbb/task_scheduler_init.h>
#include <tbb/tbb.h>

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

template <typename T>
struct count_lpdf {
  const std::vector<int>& data_;
  const T& lambda_;

  count_lpdf(const std::vector<int>& data, const T& lambda)
      : data_(data), lambda_(lambda) {}

  inline T operator()(std::size_t start, std::size_t end) const {
    const std::size_t elems = end - start;
    std::vector<int> partial_data;
    partial_data.insert(partial_data.end(), data_.begin() + start,
                        data_.begin() + end);
    return stan::math::poisson_lpmf(partial_data, lambda_);
  }
};

// without copying
template <typename T_rate>
struct fast_count_lpdf {
  typedef std::vector<int> T_n;
  const std::vector<int>& n_;
  const T_rate& lambda_;

  fast_count_lpdf(const std::vector<int>& data, const T_rate& lambda)
      : n_(data), lambda_(lambda) {}

  inline T_rate operator()(std::size_t start, std::size_t end) const {
    typedef typename stan::partials_return_type<T_n, T_rate>::type
        T_partials_return;

    T_partials_return logp(0.0);

    stan::scalar_seq_view<T_n> n_vec(n_);
    stan::scalar_seq_view<T_rate> lambda_vec(lambda_);
    // size_t size = max_size(n_, lambda_);
    std::size_t size = end - start;

    stan::math::operands_and_partials<T_rate> ops_partials(lambda_);

    for (std::size_t i = start; i != end; i++) {
      if (!(lambda_vec[i] == 0 && n_vec[i] == 0)) {
        if (stan::math::include_summand<false>::value)
          logp -= stan::math::lgamma(n_vec[i] + 1.0);
        if (stan::math::include_summand<false, T_rate>::value)
          logp += stan::math::multiply_log(n_vec[i],
                                           stan::math::value_of(lambda_vec[i]))
                  - stan::math::value_of(lambda_vec[i]);
      }

      if (!stan::is_constant_struct<T_rate>::value)
        ops_partials.edge1_.partials_[i]
            += n_vec[i] / stan::math::value_of(lambda_vec[i]) - 1.0;
    }
    return ops_partials.build(logp);
  }
};

static std::mutex cout_mutex;

struct benchmark : public ::testing::Test {
  const std::size_t elems = 10000;
  const std::size_t num_iter = 1000;
  std::vector<int> data;
  double lambda_d = 10.0;

  virtual void SetUp() {
    data.resize(elems);
    for (std::size_t i = 0; i != elems; ++i)
      data[i] = i;
  }

  virtual void TearDown() { stan::math::recover_memory(); }
};

/*
static tbb::task_scheduler_init task_scheduler(
    stan::math::internal::get_num_threads());

class ad_tape_observer : public tbb::task_scheduler_observer {
 public:
  // affinity_mask_t m_mask; // HW affinity mask to be used with an arena
  ad_tape_observer() : tbb::task_scheduler_observer() {
    std::cout << "Observer created" << std::endl;
    observe(true);  // activate the observer
  }

  ~ad_tape_observer() { std::cout << "Observer destructed" << std::endl; }

  void on_scheduler_entry(bool worker) {
    stan::math::ChainableStack::init();
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    if (worker)
      std::cout << "a worker thread is joining to do some work." << std::endl;
    else
      std::cout << "the main thread is joining to do some work." << std::endl;
  }
  void on_scheduler_exit(bool worker) {
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    if (worker)
      std::cout << "a worker thread is finished to do some work." << std::endl;
    else
      std::cout << "the main thread is finished to do some work." << std::endl;
  }
};

static ad_tape_observer tape_initializer;
*/

/*
TEST_F(benchmark, parallel_reduce_sum_d) {
  typedef boost::counting_iterator<std::size_t> count_iter;

  double lambda = lambda_d;

  count_lpdf<double> reduce_op(data, lambda);

  double poisson_lpdf = stan::math::parallel_reduce_sum(
      count_iter(0), count_iter(elems), 0.0, reduce_op);

  double poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda);

  EXPECT_FLOAT_EQ(poisson_lpdf, poisson_lpdf_ref);
}

TEST_F(benchmark, parallel_reduce_sum_v) {
  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  var lambda = lambda_d;

  count_lpdf<var> reduce_op(data, lambda);

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

TEST_F(benchmark, fast_parallel_reduce_sum_speed) {
  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  for (std::size_t i = 0; i != num_iter; ++i) {
    var lambda = lambda_d;

    fast_count_lpdf<var> reduce_op(data, lambda);

    var poisson_lpdf = stan::math::parallel_reduce_sum(
        count_iter(0), count_iter(elems), var(0.0), reduce_op);

    stan::math::grad(poisson_lpdf.vi_);
    stan::math::recover_memory();
  }
}

TEST_F(benchmark, fast_count_serial) {
  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  for (std::size_t i = 0; i != num_iter; ++i) {
    var lambda = lambda_d;

    fast_count_lpdf<var> reduce_op(data, lambda);

    var poisson_lpdf = reduce_op(0, elems);

    stan::math::grad(poisson_lpdf.vi_);
    stan::math::recover_memory();
  }
}

TEST_F(benchmark, parallel_reduce_sum_speed) {
  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  for (std::size_t i = 0; i != num_iter; ++i) {
    var lambda = lambda_d;

    count_lpdf<var> reduce_op(data, lambda);

    var poisson_lpdf = stan::math::parallel_reduce_sum(
        count_iter(0), count_iter(elems), var(0.0), reduce_op);

    stan::math::grad(poisson_lpdf.vi_);
    stan::math::recover_memory();
  }
}

TEST_F(benchmark, count_serial) {
  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  for (std::size_t i = 0; i != num_iter; ++i) {
    var lambda = lambda_d;

    count_lpdf<var> reduce_op(data, lambda);

    var poisson_lpdf = reduce_op(0, elems);

    stan::math::grad(poisson_lpdf.vi_);
    stan::math::recover_memory();
  }
}

TEST_F(benchmark, serial) {
  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  for (std::size_t i = 0; i != num_iter; ++i) {
    var lambda_ref = lambda_d;
    var poisson_lpdf_ref = stan::math::poisson_lpmf(data, lambda_ref);

    stan::math::grad(poisson_lpdf_ref.vi_);
    stan::math::recover_memory();
  }
}

TEST_F(benchmark, parallel_reduce_sum_lambda) {
  typedef boost::counting_iterator<std::size_t> count_iter;
  using stan::math::var;

  for (std::size_t i = 0; i != num_iter; ++i) {
    var lambda = lambda_d;

    var poisson_lpdf = stan::math::parallel_reduce_sum(
        count_iter(0), count_iter(elems), var(0.0),
        [&](std::size_t start, std::size_t end) {
          const std::vector<int> partial_data(data.begin() + start,
                                              data.begin() + end);
          return stan::math::poisson_lpmf(partial_data, lambda);
        });

    stan::math::grad(poisson_lpdf.vi_);
    stan::math::recover_memory();
  }
}

/*
TEST(AgradAutoDiff, parallel_for_each) {
  typedef boost::counting_iterator<int> count_iter;

  const int num_jobs = 1000;
  fun1 f;

  stan::math::var fixed_arg = 5;

  vector_v x_ref_2(num_jobs);

  for (int i = 0; i < num_jobs; ++i) {
    x_ref_2(i) = 7 + i;
  }

  auto loop_fun = [&](int i) -> vector_v {
    vector_v res(1);
    vector_v iarg(2);
    iarg << fixed_arg, x_ref_2(i);
    res(0) = f(iarg);
    return res;
  };

  vector_v parallel_result = stan::math::concatenate_row(
      stan::math::parallel_map(count_iter(0), count_iter(num_jobs), loop_fun));

  for (int i = 0; i < num_jobs; ++i) {
    vector_d x_ref(2);
    x_ref << 5, value_of(x_ref_2(i));
    double fx_ref = value_of(parallel_result(i));
    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                    fx_ref);
    vector<double> grad_fx_ref(2);
    // std::cout << "i = " << i << "; setting adjoints to zero." << std::endl;
    // stan::math::set_zero_all_adjoints_global();
    stan::math::set_zero_all_adjoints();
    // std::cout << "i = " << i << "; calling grad." << std::endl;
    stan::math::grad(parallel_result(i).vi_);
    grad_fx_ref[0] = fixed_arg.adj();
    grad_fx_ref[1] = x_ref_2(i).adj();

    EXPECT_EQ(2, grad_fx_ref.size());
    EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref[0]);
    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref[1]);
  }

  std::cout << "All good! Cleaning up..." << std::endl;
  stan::math::recover_memory();
  // stan::math::recover_memory();
}

static tbb::enumerable_thread_specific<int> arena_id;

class arena_observer : public tbb::task_scheduler_observer {
  const int id_;  // unique observer/arena id within a test
  void on_scheduler_entry(bool worker) {
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    if (worker)
      std::cout << "a worker thread is joining to do some work." << std::endl;
    else
      std::cout << "the main thread is joining to do some work." << std::endl;

    arena_id.local() = id_;

    std::cout << "this is from arena " << id_ << " thread "
              << tbb::this_task_arena::current_thread_index() << " / "
              << tbb::this_task_arena::max_concurrency() << std::endl;
  }

  void on_scheduler_exit(bool worker) {
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    if (worker)
      std::cout << "a worker thread is finished to do some work from arena "
                << id_ << ": thread "
                << tbb::this_task_arena::current_thread_index() << " / "
                << tbb::this_task_arena::max_concurrency() << std::endl;
    else
      std::cout << "the main thread is finished to do some work from arena "
                << id_ << ": thread "
                << tbb::this_task_arena::current_thread_index() << " / "
                << tbb::this_task_arena::max_concurrency() << std::endl;
  }

 public:
  arena_observer(tbb::task_arena& a, int id)
      : tbb::task_scheduler_observer(a), id_(id) {
    observe(true);
  }
};
*/
/*
TEST(AgradAutoDiff, parallel_for_each_scalar) {
  typedef boost::counting_iterator<int> count_iter;

  tbb::task_arena arena1(4);
  arena_observer obs1(arena1, 1);

  tbb::task_arena arena2(4);
  arena_observer obs2(arena2, 2);

  tbb::task& test_task = tbb::task::self();

  const int num_jobs = 1000;
  fun1 f;

  stan::math::var fixed_arg = 5;

  vector_v x_ref_2(num_jobs);

  for (int i = 0; i < num_jobs; ++i) {
    x_ref_2(i) = 7 + i;
  }

  auto loop_fun = [&](int i) -> stan::math::var {
    vector_v iarg(2);
    iarg << fixed_arg, x_ref_2(i);
    return f(iarg);
  };

  vector_v parallel_result;

  arena1.execute([&] {
    parallel_result = stan::math::concatenate_row(stan::math::parallel_map(
        count_iter(0), count_iter(num_jobs), loop_fun));
  });

  test_task.wait_for_all();

  for (int i = 0; i < num_jobs; ++i) {
    vector_d x_ref(2);
    x_ref << 5, value_of(x_ref_2(i));
    double fx_ref = value_of(parallel_result(i));
    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) * x_ref(1) + 3 * x_ref(1) * x_ref(1),
                    fx_ref);
    vector<double> grad_fx_ref(2);
    // stan::math::set_zero_all_adjoints_global();
    stan::math::set_zero_all_adjoints();
    stan::math::grad(parallel_result(i).vi_);
    grad_fx_ref[0] = fixed_arg.adj();
    grad_fx_ref[1] = x_ref_2(i).adj();

    EXPECT_EQ(2, grad_fx_ref.size());
    EXPECT_FLOAT_EQ(2 * x_ref(0) * x_ref(1), grad_fx_ref[0]);
    EXPECT_FLOAT_EQ(x_ref(0) * x_ref(0) + 3 * 2 * x_ref(1), grad_fx_ref[1]);
    // std::cout << "OK?!" << std::endl;
  }

  std::cout << "All good! Cleaning up..." << std::endl;
  stan::math::recover_memory();
  // stan::math::recover_memory();
}
 */

// test threaded AD if enabled
// more verbose old approach scheme
/*
TEST(AgradAutoDiff, parallel_for_each) {
  using stan::math::ChainableStack;
  using stan::math::var;
  typedef boost::counting_iterator<int> count_iter;
  using std::size_t;

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

  typedef ChainableStack::AutodiffStackStorage local_ad_stack_t;
  typedef std::vector<stan::math::vari*>::reverse_iterator rev_it_t;

  vector<bool> stack_is_local(num_jobs, false);
  vector<local_ad_stack_t*> stack_used(num_jobs, nullptr);
  vector<size_t> stack_starts(num_jobs);
  vector<size_t> stack_ends(num_jobs);

  std::thread::id parent_thread_id = std::this_thread::get_id();

  auto apply_f = [&](int i) -> void {
                   stack_is_local[i] = parent_thread_id ==
std::this_thread::get_id();
                   //local_ad_stack_t* local_instance = stack_is_local[i]
                   //                                   ? nullptr
                   //                                   :
&ChainableStack::instance(); local_ad_stack_t& local_instance =
ChainableStack::instance(); stack_starts[i] = local_instance.var_stack_.size();
                   fres[i] = f(x_ref_v[i]);
                   stack_ends[i] = local_instance.var_stack_.size();
                   stack_used[i] = &local_instance;
                   std::cout << "job i = " << i << "; stack size = " <<
stack_ends[i] - stack_starts[i] << std::endl;
                 };

#ifdef STAN_THREADS
  constexpr std_par::execution::parallel_unsequenced_policy exec_policy
      = std_par::execution::par_unseq;
#else
  constexpr std_par::execution::sequenced_policy exec_policy
      = std_par::execution::seq;
#endif

  std::cout << "Firing off parallel jobs" << std::endl;
  std_par::for_each(exec_policy, count_iter(0), count_iter(num_jobs), apply_f);
*/
/*
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
                //it_t begin = local_instance.var_stack_.rbegin();
                //it_t end = local_instance.var_stack_.rend();
                //for (it_t it = begin; it != end; ++it) {
                //  (*it)->chain();
                //}
              });
*/
/*
// now register all the AD tape pieces on the per-thread instances
// in the main AD tape
for (int i = 0; i < num_jobs; ++i) {
  if(!stack_is_local[i]) {
    std::cout << "non-local job i = " << i << "; stack size = " << stack_ends[i]
- stack_starts[i] << std::endl;
    stan::math::register_nested_chainablestack(*stack_used[i], stack_starts[i],
stack_ends[i]);
  }
}

for (int i = 0; i < num_jobs; ++i) {
  vector_d x_ref = x_ref_d[i];
  double fx_ref = value_of(fres[i]);
  vector<double> grad_fx_ref(2);
  //stan::math::set_zero_all_adjoints_global();
  //stan::math::set_zero_all_adjoints();
  //stan::math::grad_global(fres[i].vi_);
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

*/
