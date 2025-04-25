#ifndef STAN_TEST_UNIT_MATH_MIX_LAPLACE_UTILITY_HPP
#define STAN_TEST_UNIT_MATH_MIX_LAPLACE_UTILITY_HPP
#include <stan/math/mix.hpp>
#include <test/unit/math/mix/laplace/aki_disease_data/x1.hpp>
#include <iostream>
#include <istream>
#include <fstream>
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

/* Functions and functors used in several lgp tests. */

/////////////////////////////////////////////////////////////////////
// Covariance functions

// Function to construct spatial covariance matrix.
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> covariance(
    Eigen::Matrix<T, Eigen::Dynamic, 1> phi, int M,
    bool space_matters = false) {
  using std::pow;
  T sigma = phi[0];
  T rho = phi[1];
  double exponent;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma(M, M);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < i; j++) {
      if (space_matters) {
        exponent = i - j;
      } else {
        exponent = 1;
      }
      Sigma(i, j) = pow(rho, exponent) * sigma;
      Sigma(j, i) = Sigma(i, j);
    }
    Sigma(i, i) = sigma;
  }

  return Sigma;
}

struct spatial_covariance {
  template <typename T1, typename T2>
  Eigen::Matrix<typename stan::return_type<T1, T2>::type, Eigen::Dynamic,
                Eigen::Dynamic>
  operator()(const std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1>>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi, int M = 0) const {
    typedef typename stan::return_type<T1, T2>::type scalar;
    int space_matters = true;
    using std::pow;
    scalar sigma = phi[0];
    scalar rho = phi[1];
    double exponent;

    Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Sigma(M, M);

    for (int i = 0; i < M; i++) {
      for (int j = 0; j < i; j++) {
        if (space_matters) {
          exponent = i - j;
        } else {
          exponent = 1;
        }
        Sigma(i, j) = pow(rho, exponent) * sigma;
        Sigma(j, i) = Sigma(i, j);
      }
      Sigma(i, i) = sigma;
    }
    return Sigma;
  }
};

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
  operator()(const T1& x, std::tuple<T2, T3> arg1,
             std::ostream* msgs = nullptr) const {
    return stan::math::gp_exp_quad_cov(x, std::get<0>(arg1), std::get<1>(arg1))
           + 1e-9 * Eigen::MatrixXd::Identity(x.size(), x.size());
  }
};

// TO DO: delete this structure.
// To experiment with the prototype, provide a built-in covariance
// function. In the final version, the user will pass the covariance
// function.
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

// Naive implementation of the functor (a smarter implementation
// precomputes the covariance matrix).
struct inla_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic,
                       1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parm,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__ = 0) const {
    using stan::math::head;
    using stan::math::tail;
    using stan::math::to_vector;

    int n_groups = theta.size();
    Eigen::VectorXd n_samples = to_vector(head(dat, n_groups));
    Eigen::VectorXd sums = to_vector(tail(dat, dat.size() - n_groups));
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> Sigma
        = stan::math::test::covariance(parm, n_groups, 1);

    return sums - stan::math::elt_multiply(n_samples, stan::math::exp(theta))
           - stan::math::mdivide_left(Sigma, theta);
  }
};

// simple case where the covariance matrix is diagonal.
struct lgp_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic,
                       1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__) const {
    typedef typename stan::return_type<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> fgrad;
    int dim_theta = 2;

    Eigen::VectorXd n_samples(dim_theta);
    n_samples(0) = dat[0];
    n_samples(1) = dat[1];

    Eigen::VectorXd sums(dim_theta);
    sums(0) = dat[2];
    sums(1) = dat[3];

    return sums - stan::math::elt_multiply(n_samples, stan::math::exp(theta))
           - theta / phi(0);
  }
};

// Function to read in data for computer experiment.
// Note y and index are only required to compute the likelihood,
// although it is more efficient to do this using sufficient
// statistics.
void read_in_data(int dim_theta, int n_observations, std::string data_directory,
                  std::vector<int>& y, std::vector<int>& index,
                  std::vector<int>& sums, std::vector<int>& n_samples,
                  bool get_raw_data = false) {
  std::ifstream input_data;
  std::string dim_theta_string = std::to_string(dim_theta);
  std::string file_y = data_directory + "y_" + dim_theta_string + ".csv";
  std::string file_index
      = data_directory + "index_" + dim_theta_string + ".csv";
  std::string file_m = data_directory + "m_" + dim_theta_string + ".csv";
  std::string file_sums = data_directory + "sums_" + dim_theta_string + ".csv";

  input_data.open(file_m);
  double buffer = 0.0;
  for (int n = 0; n < dim_theta; ++n) {
    input_data >> buffer;
    n_samples[n] = buffer;
  }
  input_data.close();

  input_data.open(file_sums);
  buffer = 0.0;
  for (int n = 0; n < dim_theta; ++n) {
    input_data >> buffer;
    sums[n] = buffer;
  }
  input_data.close();

  if (get_raw_data) {
    input_data.open(file_y);
    buffer = 0.0;
    for (int n = 0; n < n_observations; ++n) {
      input_data >> buffer;
      y[n] = buffer;
    }
    input_data.close();

    input_data.open(file_index);
    buffer = 0.0;
    for (int n = 0; n < n_observations; ++n) {
      input_data >> buffer;
      index[n] = buffer;
    }
    input_data.close();
  }
}

// Overload function to read data from Aki's experiment
// using a logistic and latent Gaussian process.
void read_in_data(int dim_theta, int n_observations, std::string data_directory,
                  std::vector<double>& x1, std::vector<double>& x2,
                  std::vector<int>& y) {
  std::ifstream input_data;
  std::string file_x1 = data_directory + "x1.csv";
  std::string file_x2 = data_directory + "x2.csv";
  std::string file_y = data_directory + "y.csv";

  input_data.open(file_x1);
  double buffer = 0.0;
  for (int n = 0; n < dim_theta; ++n) {
    input_data >> buffer;
    x1[n] = buffer;
  }
  input_data.close();

  input_data.open(file_x2);
  buffer = 0.0;
  for (int n = 0; n < dim_theta; ++n) {
    input_data >> buffer;
    x2[n] = buffer;
  }
  input_data.close();

  input_data.open(file_y);
  buffer = 0.0;
  for (int n = 0; n < dim_theta; ++n) {
    input_data >> buffer;
    y[n] = buffer;
  }
  input_data.close();
}

// Overload function to read in disease mapping data.
// Same as above, but in addition include an exposure term.
void read_in_data(int dim_theta, int dim_observations,
                  std::string data_directory, std::vector<double>& x1,
                  std::vector<double>& x2, std::vector<int>& y,
                  Eigen::VectorXd& ye) {
  read_in_data(dim_theta, dim_observations, data_directory, x1, x2, y);

  std::ifstream input_data;
  std::string file_ye = data_directory + "ye.csv";

  input_data.open(file_ye);
  double buffer = 0.0;
  for (int n = 0; n < dim_theta; ++n) {
    input_data >> buffer;
    ye(n) = buffer;
  }
  input_data.close();
}

// Overload function to read in skim data.
// The covariates have a different structure.
void read_in_data(int dim_theta, int dim_observations,
                  std::string data_directory, Eigen::MatrixXd& X,
                  std::vector<int>& y, Eigen::VectorXd& lambda) {
  std::ifstream input_data;
  std::string file_y = data_directory + "y.csv";
  std::string file_X = data_directory + "X.csv";
  std::string file_lambda = data_directory + "lambda.csv";

  input_data.open(file_X);
  double buffer = 0.0;
  for (int m = 0; m < dim_theta; ++m)
    for (int n = 0; n < dim_observations; ++n) {
      input_data >> buffer;
      X(n, m) = buffer;
    }
  input_data.close();

  input_data.open(file_y);
  buffer = 0.0;
  for (int n = 0; n < dim_observations; ++n) {
    input_data >> buffer;
    y[n] = buffer;
  }
  input_data.close();

  input_data.open(file_lambda);
  buffer = 0.0;
  for (int m = 0; m < dim_theta; ++m) {
    input_data >> buffer;
    lambda[m] = buffer;
  }
}

// TODO(Steve): Should not be reading in data
// Overload function to read in gp motorcycle data.
void read_data(int dim_observations, std::string data_directory,
               std::vector<double>& x, Eigen::VectorXd& y) {
  std::ifstream input_data;
  std::string file_y = data_directory + "y_vec.csv";
  std::string file_x = data_directory + "x_vec.csv";

  input_data.open(file_y);
  double buffer = 0.0;
  y.resize(dim_observations);
  for (int i = 0; i < dim_observations; ++i) {
    input_data >> buffer;
    y(i) = buffer;
  }
  input_data.close();

  input_data.open(file_x);
  buffer = 0.0;
  x.resize(dim_observations);
  for (int i = 0; i < dim_observations; ++i) {
    input_data >> buffer;
    x[i] = buffer;
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
    data_directory = "test/unit/math/mix/laplace/aki_disease_data/";
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
    dim_phi = 2;
    phi_dbl.resize(dim_phi);
    phi_dbl << 0.3162278, 200;  // variance, length scale

    delta_lk.resize(2 * n_observations);
    y_index.resize(dim_theta);
    for (int i = 0; i < n_observations; i++) {
      delta_lk(i) = y[i];
      delta_lk(n_observations + i) = ye(i);
      y_index[i] = i;
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
  int dim_phi;
  Eigen::Matrix<double, -1, 1> phi_dbl;
  Eigen::Matrix<double, -1, 1> eta_dummy_dbl;

  Eigen::VectorXd delta_lk;
  // stan::math::poisson_log_likelihood f;
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

class laplace_count_two_dim_diag_test : public ::testing::Test {
 protected:
  void SetUp() override {
    using stan::math::algebra_solver;
    dim_theta = 2;
    y.resize(2);
    y = {1, 0};
    y_index.resize(2);
    y_index = {0, 1};

    theta_root = algebra_solver(stationary_point(), theta_0, phi, d0, di0);
    K_laplace = laplace_covariance(theta_root, phi);

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

#endif
