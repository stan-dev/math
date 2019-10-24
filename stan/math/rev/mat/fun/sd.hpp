#ifndef STAN_MATH_REV_MAT_FUN_SD_HPP
#define STAN_MATH_REV_MAT_FUN_SD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/fun/inv_sqrt.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

namespace internal {

// if x.size() = N, and x[i] = x[j] =
// then lim sd(x) -> 0 [ d/dx[n] sd(x) ] = sqrt(N) / N

inline var calc_sd(size_t size, const var* dtrs) {
  using std::sqrt;
  vari** varis = reinterpret_cast<vari**>(
      ChainableStack::instance_->memalloc_.alloc(size * sizeof(vari*)));
  double* partials = reinterpret_cast<double*>(
      ChainableStack::instance_->memalloc_.alloc(size * sizeof(double)));
  Eigen::Map<vector_vi> varis_map(varis, size);
  Eigen::Map<const vector_v> dtrs_map(dtrs, size);
  Eigen::Map<vector_d> partials_map(partials, size);

  double size_m1 = size - 1;
  varis_map = dtrs_map.vi();
  vector_d dtrs_val = dtrs_map.val();
  double mean = dtrs_val.mean();
  vector_d diff = dtrs_val.array() - mean;
  double sum_of_squares = diff.squaredNorm();
  double sd = sqrt(sum_of_squares / size_m1);

  if (sum_of_squares < 1e-20) {
    partials_map.fill(inv_sqrt(static_cast<double>(size)));
  } else {
    partials_map = diff.array() / (sd * size_m1);
  }
  return var(new stored_gradient_vari(sd, size, varis, partials));
}

}  // namespace internal

/**
 * Return the sample standard deviation of the specified standard
 * vector.  Raise domain error if size is not greater than zero.
 *
 * @param[in] v a vector
 * @return sample standard deviation of specified vector
 */
inline var sd(const std::vector<var>& v) {
  check_nonzero_size("sd", "v", v);
  if (v.size() == 1) {
    return 0;
  }
  return internal::calc_sd(v.size(), &v[0]);
}

/*
 * Return the sample standard deviation of the specified vector,
 * row vector, or matrix.  Raise domain error if size is not
 * greater than zero.
 *
 * @tparam R number of rows
 * @tparam C number of columns
 * @param[in] m input matrix
 * @return sample standard deviation of specified matrix
 */
template <int R, int C>
var sd(const Eigen::Matrix<var, R, C>& m) {
  check_nonzero_size("sd", "m", m);
  if (m.size() == 1) {
    return 0;
  }
  return internal::calc_sd(m.size(), &m(0));
}

}  // namespace math
}  // namespace stan
#endif
