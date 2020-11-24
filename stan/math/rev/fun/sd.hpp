#ifndef STAN_MATH_REV_FUN_SD_HPP
#define STAN_MATH_REV_FUN_SD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the sample standard deviation of a variable which
 * inherits from EigenBase.
 *
 * @tparam T Input type
 * @param[in] x input
 * @return sample standard deviation
 * @throw domain error  size is not greater than zero.
 */
template <typename T, require_eigen_st<is_var, T>* = nullptr>
var sd(const T& x) {
  using std::sqrt;
  using T_vi = promote_scalar_t<vari*, T>;
  using T_d = promote_scalar_t<double, T>;

  check_nonzero_size("sd", "x", x);

  if (x.size() == 1) {
    return 0.0;
  }

  vari** varis
      = ChainableStack::instance_->memalloc_.alloc_array<vari*>(x.size());
  double* partials
      = ChainableStack::instance_->memalloc_.alloc_array<double>(x.size());
  Eigen::Map<T_vi> varis_map(varis, x.rows(), x.cols());
  Eigen::Map<T_d> partials_map(partials, x.rows(), x.cols());

  varis_map = x.vi();
  T_d dtrs_val = x.val();
  double mean = dtrs_val.mean();
  T_d diff = dtrs_val.array() - mean;
  double sum_of_squares = diff.squaredNorm();
  double size_m1 = x.size() - 1;
  double sd = sqrt(sum_of_squares / size_m1);

  if (sum_of_squares < 1e-20) {
    partials_map.fill(inv_sqrt(static_cast<double>(x.size())));
  } else {
    partials_map = diff.array() / (sd * size_m1);
  }
  return var(new stored_gradient_vari(sd, x.size(), varis, partials));
}

/**
 * Return the sample standard deviation of the var_value matrix
 *
 * @tparam T Input type
 * @param[in] x input matrix
 * @return sample standard deviation of specified matrix
 * @throw domain error  size is not greater than zero.
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
var sd(const T& x) {
  check_nonzero_size("sd", "x", x);

  if (x.size() == 1) {
    return 0.0;
  }

  auto arena_diff = to_arena((x.val().array() - x.val().mean()).matrix());
  double sum_of_squares = arena_diff.squaredNorm();
  double sd = std::sqrt(sum_of_squares / (x.size() - 1));

  return make_callback_vari(sd, [x, arena_diff](const auto& res) mutable {
    x.adj() += (res.adj() / (res.val() * (x.size() - 1))) * arena_diff;
  });
}

/**
 * Return the sample standard deviation of the specified std vector, column
 * vector, row vector, matrix, or std vector of any of these types.
 *
 * @tparam T Input type
 * @param[in] m input matrix
 * @return sample standard deviation
 * @throw domain error  size is not greater than zero.
 */
template <typename T, require_std_vector_st<is_var, T>* = nullptr>
auto sd(const T& m) {
  return apply_vector_unary<T>::reduce(m, [](const auto& x) { return sd(x); });
}

}  // namespace math
}  // namespace stan
#endif
