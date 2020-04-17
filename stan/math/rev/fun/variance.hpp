#ifndef STAN_MATH_REV_FUN_VARIANCE_HPP
#define STAN_MATH_REV_FUN_VARIANCE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

inline var calc_variance(size_t size, const var* dtrs) {
  vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(size);
  double* partials
      = ChainableStack::instance_->memalloc_.alloc_array<double>(size);

  Eigen::Map<const vector_v> dtrs_map(dtrs, size);
  Eigen::Map<vector_vi>(varis, size) = dtrs_map.vi();
  vector_d dtrs_vals = dtrs_map.val();

  vector_d diff = dtrs_vals.array() - dtrs_vals.mean();
  double size_m1 = size - 1;
  Eigen::Map<vector_d>(partials, size) = 2 * diff.array() / size_m1;
  double variance = diff.squaredNorm() / size_m1;

  return var(new stored_gradient_vari(variance, size, varis, partials));
}

}  // namespace internal

/**
 * Return the sample variance of the specified standard
 * vector.  Raise domain error if size is not greater than zero.
 *
 * @param[in] v a vector
 * @return sample variance of specified vector
 */
 template <typename StdVec, require_std_vector_t<StdVec>* = nullptr,
   require_vt_var<StdVec>* = nullptr>
inline auto variance(StdVec&& v) {
  check_nonzero_size("variance", "v", v);
  if (v.size() == 1) {
    return var{0};
  }
  return var{internal::calc_variance(v.size(), &v[0])};
}

/**
 * Return the sample variance of the specified vector, row vector,
 * or matrix.  Raise domain error if size is not greater than
 * zero.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param[in] m input matrix
 * @return sample variance of specified matrix
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
  require_vt_var<EigMat>* = nullptr>
auto variance(EigMat&& m) {
  using ref_inner = const typename std::decay_t<EigMat>::PlainObject;
  check_nonzero_size("variance", "m", m);
  if (m.size() == 1) {
    return var{0};
  }

  const Eigen::Ref<ref_inner, Eigen::Aligned16, Eigen::Stride<0,0>>& mat = m;
  return var{internal::calc_variance(m.size(), mat.data())};
}

}  // namespace math
}  // namespace stan
#endif
