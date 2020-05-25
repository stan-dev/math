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
 * Return the sample standard deviation of the specified std vector, column
 * vector, row vector, or matrix.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param[in] m input matrix
 * @return sample standard deviation of specified matrix
 * @throw domain error  size is not greater than zero.
 */
template <typename T, require_container_vt<is_var, T>* = nullptr>
var sd(const T& m) {
  check_nonzero_size("sd", "m", m);
  if (m.size() == 1) {
    return 0;
  }

  return apply_vector_unary<T>::reduce(m, [](const auto& dtrs_map) {
    using std::sqrt;
    using T_map = std::decay_t<decltype(dtrs_map)>;
    using T_vi = promote_scalar_t<vari*, T_map>;
    using T_d = promote_scalar_t<double, T_map>;
    vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
        dtrs_map.size());
    double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
        dtrs_map.size());
    Eigen::Map<T_vi> varis_map(varis, dtrs_map.rows(), dtrs_map.cols());
    Eigen::Map<T_d> partials_map(partials, dtrs_map.rows(), dtrs_map.cols());

    varis_map = dtrs_map.vi();
    T_d dtrs_val = dtrs_map.val();
    double mean = dtrs_val.mean();
    T_d diff = dtrs_val.array() - mean;
    double sum_of_squares = diff.squaredNorm();
    double size_m1 = dtrs_map.size() - 1;
    double sd = sqrt(sum_of_squares / size_m1);

    if (sum_of_squares < 1e-20) {
      partials_map.fill(inv_sqrt(static_cast<double>(dtrs_map.size())));
    } else {
      partials_map = diff.array() / (sd * size_m1);
    }
    return var(new stored_gradient_vari(sd, dtrs_map.size(), varis, partials));
  });
}

}  // namespace math
}  // namespace stan
#endif
