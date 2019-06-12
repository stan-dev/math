#include <gtest/gtest.h>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

/* Functors and functions for the performance tests. */

// evolution operator for a one compartment model
// Signature needs to be the same as the one used for the function.
template <typename T0, typename T1>
inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                     Eigen::Dynamic, 1>
oneCpt_evolution(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& y,
                 const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
                 double dt,
                 int n_patients) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  typedef typename stan::return_type<T0, T1>::type scalar;

  Matrix<scalar, Dynamic, 1> state(2 * n_patients);

  // extract population parameters
  T1 k1_pop = theta(0);
  T1 k2_pop = theta(1);

  // compute states for all patients
  for (int i = 0; i < n_patients; i++) {
    T1 k1 = k1_pop * theta(2*i + 2);
    T1 k2 = k2_pop * theta(2*i + 3);

    state(2*i) = y(2*i) * exp(-k1 * dt);
    state(2*i + 1) = y(2*i) * k1 / (k2 - k1) * (exp(-k1 * dt)
                       - exp(-k2 * dt)) + y(2*i + 1) * exp(-k2 * dt);
  }

  return state;
}

// functor to encode system of algebraic equations
struct oneCpt_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator() (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& y,
           const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int,
           std::ostream* pstream__ = 0) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;

    double dose = dat[0];
    double dt = dat[1];
    int n_patients = dat_int[0];

    Matrix<T0, Dynamic, 1> y_init = y;

    for (int i = 0; i < n_patients; i++)
      y_init(2*i) += dose; 

    return oneCpt_evolution(y_init, theta, dt, n_patients) - y;
  }
};
