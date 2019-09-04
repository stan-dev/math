#ifndef STAN_MATH_REV_MAT_FUNCTOR_CHECK_FLAG_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_CHECK_FLAG_HPP

namespace stan {
namespace math {

// TO DO (charlesm93): use this function inside cvodes
// integrator.
/** 
 * Throws an exception when a Sundial function fails 
 * (i.e. returns a negative flag)
 */
inline void check_flag(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func_name << " failed with error flag " << flag << ".";
    throw std::runtime_error(ss.str());
  }
}

/**
 * Throws an exception message when the function KINSol()
 * (call to the solver) fails. When the exception is caused
 * by a tuning parameter the user controls, gives a specific
 * error.
 */
inline void check_flag_kinsol(int flag, long int max_num_steps) {
  std::ostringstream ss;
  if (flag == -6) {
    ss << "algebra_solver: max number of iterations: "
       << max_num_steps
       << " exceeded.";
    throw boost::math::evaluation_error(ss.str());
  } else if (flag < 0) {
    ss << "algebra_solver failed with error flag " << flag << ".";
    throw boost::math::evaluation_error(ss.str());
  }
}

}  // namespace math
}  // namespace stan
#endif
