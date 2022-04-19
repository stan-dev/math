#ifndef STAN_MATH_REV_FUNCTOR_ALGEBRA_SYSTEM_HPP
#define STAN_MATH_REV_FUNCTOR_ALGEBRA_SYSTEM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/jacobian.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * A structure which gets passed to Eigen's dogleg
 * algebraic solver.
 *
 * @tparam T scalar type of independent variable.
 * @tparam NX number of rows
 * @tparam NY number of columns
 */
template <typename T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct nlo_functor {
  const int m_inputs, m_values;

  nlo_functor() : m_inputs(NX), m_values(NY) {}

  nlo_functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

/**
 * A functor with the required operators to call Eigen's algebraic solver.
 *
 * @tparam S wrapper around the algebraic system functor. Has the
 * signature required for jacobian (i.e takes only one argument).
 * @tparam F algebraic system functor
 * @tparam T0 scalar type for unknowns
 * @tparam T1 scalar type for auxiliary parameters
 */
template <typename S>
struct hybrj_functor_solver : nlo_functor<double> {
  /** Wrapper around algebraic system */
  S fs_;

  /** Jacobian of algebraic function wrt unknowns */
  Eigen::MatrixXd J_;

  explicit hybrj_functor_solver(const S& fs) : fs_(fs) {}

  /**
   * Computes the value the algebraic function, f, when pluging in the
   * independent variables, and the Jacobian w.r.t unknowns.
   *
   * @param [in] iv independent variables
   * @param [in, out] fvec value of algebraic function when plugging in iv.
   */
  int operator()(const Eigen::VectorXd& iv, Eigen::VectorXd& fvec) {
    jacobian(fs_, iv, fvec, J_);
    return 0;
  }

  /**
   * Assign the Jacobian to fjac.
   *
   * @param [in] iv independent variables.
   * @param [in, out] fjac matrix container for jacobian
   */
  int df(const Eigen::VectorXd& iv, Eigen::MatrixXd& fjac) const {
    fjac = J_;
    return 0;
  }

  /**
   * Performs the same task as the operator(), but returns the
   * Jacobian, instead of saving it inside an argument
   * passed by reference.
   * @param [in] iv independent variable.
   */
  Eigen::MatrixXd get_jacobian(const Eigen::VectorXd& iv) {
    Eigen::VectorXd fvec;
    jacobian(fs_, iv, fvec, J_);
    return J_;
  }

  /**
   * Performs the same task as df(), but returns the value of
   * algebraic function, instead of saving it inside an
   * argument passed by reference.
   * @tparam [in] iv independent variable.
   */
  Eigen::VectorXd get_value(const Eigen::VectorXd& iv) const { return fs_(iv); }
};

// TODO(jgaeb): Remove this when the chain method of the fixed point solver is
// updated.
template <typename T1, typename T2>
void algebra_solver_check(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                          const Eigen::Matrix<T2, Eigen::Dynamic, 1> y,
                          const std::vector<double>& dat,
                          const std::vector<int>& dat_int,
                          double function_tolerance,
                          long int max_num_steps) {  // NOLINT(runtime/int)
  check_nonzero_size("algebra_solver", "initial guess", x);
  check_finite("algebra_solver", "initial guess", x);
  check_finite("algebra_solver", "parameter vector", y);
  check_finite("algebra_solver", "continuous data", dat);
  check_finite("algebra_solver", "integer data", dat_int);

  check_nonnegative("algebra_solver", "function_tolerance", function_tolerance);
  check_positive("algebra_solver", "max_num_steps", max_num_steps);
}

}  // namespace math
}  // namespace stan

#endif
