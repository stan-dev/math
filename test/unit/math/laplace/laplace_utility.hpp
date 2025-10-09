#ifndef STAN_TEST_UNIT_MATH_MIX_LAPLACE_UTILITY_HPP
#define STAN_TEST_UNIT_MATH_MIX_LAPLACE_UTILITY_HPP
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/aki_disease_data/x1.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <gtest/gtest.h>

namespace stan {
namespace math {
namespace test {

struct laplace_issue {
  int solver_num;
  int max_steps_line_search;
  int hessian_block_size;
  constexpr laplace_issue(int solv, int max_steps, int hess_block)
      : solver_num(solv),
        max_steps_line_search(max_steps),
        hessian_block_size(hess_block) {}
  constexpr bool operator==(const laplace_issue& other) const {
    return solver_num == other.solver_num
           && max_steps_line_search == other.max_steps_line_search
           && hessian_block_size == other.hessian_block_size;
  }
};

namespace internal {
template <typename T>
struct is_pair : std::false_type {};

template <typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type {};
}  // namespace internal

template <typename T>
static constexpr bool is_pair_v = internal::is_pair<std::decay_t<T>>::value;

enum class LaplaceFailures {
  HessianFailure = 0,
  SqrtDNE = 1,
  NaNTheta = 2,
  IterExceeded = 3,
  Other = 4,
  None = 5
};
inline std::string to_string(LaplaceFailures value) {
  switch (value) {
    case LaplaceFailures::HessianFailure:
      return "LaplaceFailures::HessianFailure";
    case LaplaceFailures::SqrtDNE:
      return "LaplaceFailures::SqrtDNE";
    case LaplaceFailures::NaNTheta:
      return "LaplaceFailures::NaNTheta";
    case LaplaceFailures::IterExceeded:
      return "LaplaceFailures::IterExceeded";
    case LaplaceFailures::None:
      return "LaplaceFailures::None";
    default:
      return "LaplaceFailures::Other";
  }
}

template <typename T>
inline auto err_to_laplace_failure(T&& e) {
  if (std::string(e.what()).find("positive") != std::string::npos) {
    return LaplaceFailures::HessianFailure;
  } else if (std::string(e.what()).find("schur") != std::string::npos) {
    return LaplaceFailures::SqrtDNE;
  } else if (std::string(e.what()).find("NaN") != std::string::npos) {
    return LaplaceFailures::NaNTheta;
  } else if (std::string(e.what()).find("iteration") != std::string::npos) {
    return LaplaceFailures::IterExceeded;
  } else {
    return LaplaceFailures::Other;
  }
  return LaplaceFailures::None;
}

template <typename T1, typename T2>
inline constexpr auto flag_test(T1&& known_issues, T2&& test_params) {
  if constexpr (is_pair_v<decltype(known_issues[0])>) {
    for (auto&& issue : known_issues) {
      if (issue.first == test_params) {
        return issue.second;
      }
    }
    return LaplaceFailures::None;
  } else {
    for (auto&& issue : known_issues) {
      if (issue == test_params) {
        return true;
      }
    }
    return false;
  }
}

template <typename T1>
inline constexpr auto flag_test(T1&& known_issues, int solver_num,
                                int max_steps_line_search,
                                int hessian_block_size) {
  return flag_test(
      known_issues,
      laplace_issue{solver_num, max_steps_line_search, hessian_block_size});
}

struct squared_kernel_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> operator()(
      const T2& x, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
      const std::vector<double>& delta, const std::vector<int>& delta_int,
      std::ostream* msgs = nullptr) const {
    return stan::math::gp_exp_quad_cov(x, phi(0), phi(1))
           + 1e-9 * Eigen::MatrixXd::Identity(x.size(), x.size());
  }
  template <typename T1, typename T2, typename T3>
  Eigen::Matrix<return_type_t<T1, T2, T3>, Eigen::Dynamic, Eigen::Dynamic>
  operator()(const T1& x, const T2& arg1, const T3& arg2,
             std::ostream* msgs = nullptr) const {
    return stan::math::gp_exp_quad_cov(x, arg1, arg2)
           + 1e-9 * Eigen::MatrixXd::Identity(x.size(), x.size());
  }
  template <typename T1, typename T2, typename T3>
  Eigen::Matrix<return_type_t<T1, T2, T3>, Eigen::Dynamic, Eigen::Dynamic>
  operator()(const T1& x, const std::tuple<T2, T3>& arg1,
             std::ostream* msgs = nullptr) const {
    return stan::math::gp_exp_quad_cov(x, std::get<0>(arg1), std::get<1>(arg1))
           + 1e-9 * Eigen::MatrixXd::Identity(x.size(), x.size());
  }
  template <typename T1, typename T2, typename T3>
  Eigen::Matrix<return_type_t<T1, T2, T3>, Eigen::Dynamic, Eigen::Dynamic>
  operator()(const T1& x, const std::vector<std::tuple<T2, T3>>& arg1,
             std::ostream* msgs = nullptr) const {
    return stan::math::gp_exp_quad_cov(x, std::get<0>(arg1[0]),
                                       std::get<1>(arg1[0]))
           + 1e-9 * Eigen::MatrixXd::Identity(x.size(), x.size());
  }
};

struct sqr_exp_kernel_functor {
  template <typename T1, typename T2, typename T3>
  auto operator()(const T1& x, const T2& alpha, const T3& rho,
                  std::ostream* msgs = nullptr) const {
    constexpr double jitter = 1e-8;
    Eigen::Matrix<return_type_t<T1, T2, T3>, Eigen::Dynamic, Eigen::Dynamic>
        kernel = stan::math::gp_exp_quad_cov(x, alpha, rho);
    for (int i = 0; i < kernel.cols(); i++)
      kernel(i, i) += jitter;

    return kernel;
  }
};

struct stationary_point {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic,
                       1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parms,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__ = 0) const {
    Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic, 1>
        z(2);
    z(0) = 1 - exp(theta(0)) - theta(0) / (parms(0) * parms(0));
    z(1) = 0 - exp(theta(1)) - theta(1) / (parms(1) * parms(1));
    return z;
  }
};

struct diagonal_kernel_functor {
  template <typename T1, typename T2>
  auto operator()(const T1& alpha, const T2& rho,
                  std::ostream* msgs = nullptr) const {
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
    K(0, 0) = alpha * alpha;
    K(1, 1) = rho * rho;
    K(0, 1) = 0;
    K(1, 0) = 0;
    return K;
  }
};

template <typename F, typename ThetaVec>
inline void run_solver_grid(F&& body, ThetaVec&& theta_0) {
  constexpr std::array solver_nums{1, 2, 3};            // [1, 3]
  constexpr std::array hessian_block_sizes{1, 2, 3};    // [1, 2]
  constexpr std::array max_steps_line_searches{0, 10};  // 0, 10
  for (int solver : solver_nums) {
    for (int hblock : hessian_block_sizes) {
      for (int ls_steps : max_steps_line_searches) {
        if (theta_0.size() % hblock != 0) {
          std::cerr << "[          ] [ INFO ]"
                    << " Skipping test for hessian of size " << theta_0.size()
                    << " with hessian block size of " << hblock << std::endl;
          continue;
        }
        try {
          std::forward<F>(body)(solver, hblock, ls_steps, theta_0);
        } catch (const std::exception& e) {
          ADD_FAILURE() << "Exception: " << e.what();
        }
        if (::testing::Test::HasFailure()) {
          std::cout << "----------" << std::endl;
          std::cout << "solver_num: " << solver << std::endl;
          std::cout << "hessian_block_size: " << hblock << std::endl;
          std::cout << "max_steps_line_search: " << ls_steps << std::endl;
        }
      }
    }
  }
}

template <typename T1, typename T2>
Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> laplace_covariance(
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_root,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& phi) {
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
  K(0, 0) = 1 / (stan::math::exp(theta_root(0)) + 1 / (phi(0) * phi(0)));
  K(1, 1) = 1 / (stan::math::exp(theta_root(1)) + 1 / (phi(1) * phi(1)));
  K(0, 1) = 0;
  K(1, 0) = 0;
  return K;
}

/**
 * Helper function for printing out adjoints
 */
template <typename Output, require_t<is_any_var_scalar<Output>>* = nullptr>
inline void print_adjoint(Output&& output) {
  if constexpr (is_tuple_v<Output>) {
    std::cout << "tuple adj\n";
    return stan::math::for_each(
        [](auto&& output_i) { return print_adjoint(output_i); }, output);
  } else if constexpr (is_std_vector_v<Output>) {
    if constexpr (is_var_v<value_type_t<Output>>) {
      Eigen::Map<const Eigen::Matrix<var, -1, -1>> map_x(output.data(),
                                                         output.size());
      std::cout << "eigen adj: \n" << map_x.adj() << std::endl;
    } else {
      std::cout << "stdvec adjoint\n";
      for (int i = 0; i < output.size(); ++i) {
        print_adjoint(output[i]);
      }
    }
  } else if constexpr (is_eigen_v<Output>) {
    std::cout << "adj: \n" << output.adj() << std::endl;
  } else if constexpr (is_stan_scalar_v<Output>) {
    std::cout << "adj: " << output.adj() << std::endl;
  } else {
    static_assert(sizeof(Output*) == 0,
                  "INTERNAL ERROR:(laplace_marginal_lpdf) print_adjoint was "
                  "not able to deduce the actiopns needed for the given type. "
                  "This is an internal error, please report it: "
                  "https://github.com/stan-dev/math/issues");
  }
}

}  // namespace test
}  // namespace math
}  // namespace stan

//////////////////////////////////////////////////////////////////////////

class laplace_disease_map_test : public ::testing::Test {
  // Based on (Vanhatalo, Pietilainen and Vethari, 2010). See
  // https://research.cs.aalto.fi/pml/software/gpstuff/demo_spatial1.shtml
 protected:
  void SetUp() override {
    dim_theta = 911;
    n_observations = 911;
    x1 = stan::test::laplace::disease::x1;
    x2 = stan::test::laplace::disease::x2;
    y = stan::test::laplace::disease::y;
    ye = stan::test::laplace::disease::ye;

    dim_x = 2;
    x.resize(dim_theta);
    for (int i = 0; i < dim_theta; i++) {
      Eigen::VectorXd coordinate(dim_x);
      coordinate << x1[i], x2[i];
      x[i] = coordinate;
    }

    // one observation per group
    n_samples.resize(dim_theta);
    for (int i = 0; i < dim_theta; i++)
      n_samples[i] = 1;

    theta_0 = Eigen::VectorXd::Zero(dim_theta);
    mean = Eigen::VectorXd::Zero(dim_theta);
    dim_phi = 2;
    phi_dbl.resize(dim_phi);
    phi_dbl << 0.3162278, 200;  // variance, length scale

    delta_lk.resize(2 * n_observations);
    y_index.resize(dim_theta);
    for (int i = 0; i < n_observations; i++) {
      delta_lk(i) = y[i];
      delta_lk(n_observations + i) = ye(i);
      y_index[i] = i + 1;
    }
  }

  int dim_theta;
  int n_observations;
  std::string data_directory;
  std::vector<double> x1, x2;
  std::vector<int> y;
  Eigen::VectorXd ye;
  int dim_x;
  std::vector<Eigen::VectorXd> x;
  std::vector<int> y_index;
  std::vector<int> n_samples;
  std::vector<double> delta;
  std::vector<int> delta_int;

  Eigen::VectorXd theta_0;
  Eigen::VectorXd mean;
  int dim_phi;
  Eigen::Matrix<double, -1, 1> phi_dbl;
  Eigen::Matrix<double, -1, 1> eta_dummy_dbl;

  Eigen::VectorXd delta_lk;
  // stan::math::poisson_log_likelihood f;
};

class laplace_count_two_dim_diag_test : public ::testing::Test {
 protected:
  void SetUp() override {
    using stan::math::algebra_solver;
    dim_theta = 2;
    y.resize(2);
    y = {1, 0};
    y_index.resize(2);
    y_index = {1, 2};

    theta_root = algebra_solver(stan::math::test::stationary_point(), theta_0,
                                phi, d0, di0);
    K_laplace = stan::math::test::laplace_covariance(theta_root, phi);

    rng.seed(1954);
    theta_benchmark = stan::math::multi_normal_rng(theta_root, K_laplace, rng);

    tol = 1e-3;
    n_sim = 5e5;
  }

  int dim_theta;
  Eigen::VectorXd theta_0{{1, 1}};
  Eigen::VectorXd theta_root;
  Eigen::VectorXd phi{{3, 2}};
  std::vector<int> y;
  std::vector<int> y_index;
  Eigen::VectorXd ye{{1, 1}};
  std::vector<double> d0;
  std::vector<int> di0;
  Eigen::MatrixXd K_laplace;
  Eigen::MatrixXd theta_benchmark;
  boost::random::mt19937 rng;
  double tol;
  int n_sim;
};
#ifdef DEBUG_LAPLACE
static bool write_init_json = true;
static int err_iter = 0;

// Custom event listener that logs test failures
class LoggingTestListener : public ::testing::EmptyTestEventListener {
 public:
  std::string current_test_name_;
  int solver_num{0};
  int hessian_block_size{0};
  int max_steps_line_search{0};

  // Called after an assertion results in a failure.
  void OnTestPartResult(const ::testing::TestPartResult& result) override {
    if (result.failed()) {
      std::ofstream ofs;
      // On first failure, open file in truncation mode and write header
      if (write_init_json) {
        ofs.open("failure_log.json", std::ios::out);
        ofs << "{\"error\": {\n";
        write_init_json = false;
      } else {
        // For subsequent failures, open in append mode and add a comma
        // separator
        ofs.open("failure_log.json", std::ios::app);
        ofs << ", \n";
      }
      ofs << "\"" << err_iter << "\": {";
      err_iter++;
      std::string result_summary = result.summary();
      boost::replace_all(result_summary, "\"", "\\\"");
      boost::replace_all(result_summary, "\n", "\\n");
      // Retrieve the current test information.
      const ::testing::TestInfo* test_info
          = ::testing::UnitTest::GetInstance()->current_test_info();
      std::string test_name;
      if (test_info) {
        // For Google Test 1.10.0 or later
        test_name = std::string(test_info->test_suite_name()) + "."
                    + test_info->name();
        // For older versions, use:
        // test_name = std::string(test_info->test_case_name()) + "." +
        // test_info->name();
      }
      ofs << "\"test\": \"" << test_name << "\", ";
      ofs << "\"solver_num\": " << solver_num << ", ";
      ofs << "\"hessian_block_size\": " << hessian_block_size << ", ";
      ofs << "\"max_steps_line_search\": " << max_steps_line_search << ", ";
      ofs << "\"failure\": \"" << result_summary << "\"}";
    }
  }

  // Called after all tests have ended.
  void OnTestProgramEnd(const ::testing::UnitTest& /*unit_test*/) override {
    if (!write_init_json) {  // Only if at least one failure was logged
      std::ofstream ofs("failure_log.json", std::ios::app);
      ofs << "}}";  // Close the JSON object
    }
  }
};

class laplace_test_listen : public ::testing::Test {
 public:
  virtual void AllowSetup() {}
  bool setup_once{true};
  LoggingTestListener* logger{new LoggingTestListener{}};

 protected:
  virtual void AllowSetup() {
    if (setup_once) {
      ::testing::TestEventListeners& listeners
          = ::testing::UnitTest::GetInstance()->listeners();
      listeners.Append(logger);
      setup_once = false;
    }
  }
  void SetUp() override { this->AllowSetup(); }
  virtual ~laplace_test_listen() {
    ::testing::TestEventListeners& listeners
        = ::testing::UnitTest::GetInstance()->listeners();
    listeners.Release(logger);
    delete logger;
  }
}
#endif

#endif
