#ifndef STAN_MATH_REV_CORE_VAR_VALUE_FWD_DECLARE_HPP
#define STAN_MATH_REV_CORE_VAR_VALUE_FWD_DECLARE_HPP

namespace stan {
namespace math {
// forward declaration of var
template <typename T, typename = void>
class var_value;

/**
 * Equivalent to `Eigen::Matrix`, except that the data is stored on AD stack.
 * That makes these objects trivially destructible and usable in `vari`s.
 *
 * @tparam MatrixType Eigen matrix type this works as (`MatrixXd`, `VectorXd`,
 * ...)
 */
template <typename MatrixType, typename = void>
class arena_matrix;
}  // namespace math
}  // namespace stan
#endif
