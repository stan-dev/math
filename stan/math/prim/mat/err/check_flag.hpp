#ifndef STAN_MATH_REV_MAT_FUNCTOR_CHECK_FLAG_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_CHECK_FLAG_HPP

namespace stan {
namespace math {

// TO DO (charlesm93): use this function inside cvodes
// integrator.
inline void check_flag(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func_name << " failed with error flag " << flag;
    throw std::runtime_error(ss.str());
  }
}

}  // namespace math
}  // namespace stan
#endif
