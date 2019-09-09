#include <iostream>
#include <istream>
#include <fstream>
 
/* Functions and functors used in several lgp tests. */

/////////////////////////////////////////////////////////////////////
// Covariance functions

// Function to construct spatial covariance matrix.
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
covariance (Eigen::Matrix<T, Eigen::Dynamic, 1> phi, int M,
            bool space_matters = false) {
  using std::pow;
  T sigma = phi[0];
  T rho = phi[1];
  double exponent;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma(M, M);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < i; j++) {
      if (space_matters) {exponent = i - j;} else {exponent = 1;}
      Sigma(i, j) = pow(rho, exponent) * sigma;
      Sigma(j, i) = Sigma(i, j);
    }
    Sigma(i, i) = sigma;
  }

  return Sigma;
}

struct spatial_covariance {
  template <typename T1, typename T2>
  Eigen::Matrix<typename stan::return_type<T1, T2>::type,
                Eigen::Dynamic, Eigen::Dynamic>
  operator() (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
           const std::vector<Eigen::Matrix<T2, 
                                           Eigen::Dynamic, 1>>& x,
                                           int M = 0) const {
    typedef typename stan::return_type<T1, T2>::type scalar;
    int space_matters = true;
    using std::pow;
    scalar sigma = phi[0];
    scalar rho = phi[1];
    double exponent;
    
    Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Sigma(M, M);
    
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < i; j++) {
        if (space_matters) {exponent = i - j;} else {exponent = 1;}
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
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
  operator() (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
           const T2& x, int M = 0) const {
    return stan::math::gp_exp_quad_cov(x, phi(0), phi(1))
    + 1e-9 * Eigen::MatrixXd::Identity(x.size(), x.size());
  }
};

// Naive implementation of the functor (a smarter implementation
// precomputes the covariance matrix).
struct inla_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator() (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
           const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parm,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int,
           std::ostream* pstream__ = 0) const {
    using stan::math::to_vector;
    using stan::math::head;
    using stan::math::tail;

    int n_groups = theta.size();
    Eigen::VectorXd n_samples = to_vector(head(dat, n_groups));
    Eigen::VectorXd sums = to_vector(tail(dat, dat.size() - n_groups));
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
      Sigma = covariance(parm, n_groups, 1);

    return sums - stan::math::elt_multiply(n_samples, stan::math::exp(theta)) -
      stan::math::mdivide_left(Sigma, theta);
  }
};

// simple case where the covariance matrix is diagonal.
struct lgp_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator ()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
            const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
            const std::vector<double>& dat,
            const std::vector<int>& dat_int,
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

    return sums - stan::math::elt_multiply(n_samples,
                                           stan::math::exp(theta))
      - theta / phi(0);
  }
};

// Function to read in data for computer experiment.
// Note y and index are only required to compute the likelihood,
// although it is more efficient to do this using sufficient
// statistics.
void read_in_data (int dim_theta,
                   int n_observations,
                   std::string data_directory,
                   std::vector<int>& y,
                   std::vector<int>& index,
                   std::vector<int>& sums,
                   std::vector<int>& n_samples,
                   bool get_raw_data = false) {
  std::ifstream input_data;
  std::string dim_theta_string = std::to_string(dim_theta);
  std::string file_y = data_directory + "y_" + dim_theta_string + ".csv";
  std::string file_index = data_directory + "index_" +
    dim_theta_string + ".csv";
  std::string file_m = data_directory + "m_" + dim_theta_string + ".csv";
  std::string file_sums = data_directory + "sums_" +
    dim_theta_string + ".csv";

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
void read_in_data (int dim_theta,
                   int n_observations,
                   std::string data_directory,
                   std::vector<double>& x1,
                   std::vector<double>& x2,
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
void read_in_data(int dim_theta,
                  int dim_observations,
                  std::string data_directory,
                  std::vector<double>& x1,
                  std::vector<double>& x2,
                  std::vector<int>& y,
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
