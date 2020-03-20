#ifndef STAN_MATH_REV_CORE_SAVE_VARIS_HPP
#define STAN_MATH_REV_CORE_SAVE_VARIS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

#include <vector>

namespace stan {
namespace math {

/**
 * Store a reference to `vari` used in `reduce_sum`.
 * @tparam Pargs Types to be forwarded to recursive call of `save_varis`.
 * @param[in, out] dest Pointer to set of varis.
 * @param[in] x A `var` whose `vari` will be stores
 * @param[in] args Further arguments forwarded to recursive call.
 */
template <typename... Pargs>
inline vari** save_varis(vari** dest, const var& x, Pargs&&... args) {
  *dest = x.vi_;
  return save_varis(dest + 1, std::forward<Pargs>(args)...);
}

/**
 * Store a reference to an std vector's `vari` used in `reduce_sum`.
 * @tparam VarVec Type of standard vector holding vars.
 * @tparam Pargs Types to be forwarded to recursive call of `save_varis`.
 * @param[in, out] dest Pointer to set of varis.
 * @param[in] x A standard vector whose vari's are stored.
 * @param[in] args Further arguments forwarded to recursive call.
 */
template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr,
          typename... Pargs>
inline vari** save_varis(vari** dest, VarVec&& x, Pargs&&... args) {
  using write_map = Eigen::Map<Eigen::Matrix<vari*, -1, 1>>;
  using read_map = Eigen::Map<const Eigen::Matrix<var, -1, 1>>;
  write_map(dest, x.size(), 1) = read_map(x.data(), x.size(), 1).vi();
  return save_varis(dest + x.size(), std::forward<Pargs>(args)...);
}

/**
 * Store a reference to an std vector of containers's `vari` used in
 * `reduce_sum`.
 *
 * @tparam VecContainer Type of standard vector holding containers of vars.
 * @tparam Pargs Types to be forwarded to recursive call of `save_varis`.
 * @param[in, out] dest Pointer to set of varis.
 * @param[in] x A standard vector whose container's varis are stored.
 * @param[in] args Further arguments forwarded to recursive call.
 */
template <typename VecContainer,
          require_std_vector_st<is_var, VecContainer>* = nullptr,
          require_std_vector_vt<is_container, VecContainer>* = nullptr,
          typename... Pargs>
inline vari** save_varis(vari** dest, VecContainer&& x, Pargs&&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest = save_varis(dest, x[i]);
  }
  return save_varis(dest, std::forward<Pargs>(args)...);
}

/**
 * Store a reference to an Eigen types `vari` used in `reduce_sum`.
 *
 * @tparam VarVec Type derived from `EigenBase` holding vars.
 * @tparam Pargs Types to be forwarded to recursive call of `save_varis`.
 * @param[in, out] dest Pointer to set of varis.
 * @param[in] x An `EigenBase` whose vari's are stored.
 * @param[in] args Further arguments forwarded to recursive call.
 */
template <typename EigT, require_eigen_vt<is_var, EigT>* = nullptr,
          typename... Pargs>
inline vari** save_varis(vari** dest, EigT&& x, Pargs&&... args) {
  using mat_t = std::decay_t<EigT>;
  using write_map = Eigen::Map<
      Eigen::Matrix<vari*, mat_t::RowsAtCompileTime, mat_t::ColsAtCompileTime>>;
  write_map(dest, x.rows(), x.cols()) = x.vi();
  return save_varis(dest + x.size(), std::forward<Pargs>(args)...);
}

/**
 * Specialization to skip over arithmetic types and containers of arithmetic
 * types.
 *
 * @tparam Arith either an arithmetic type or container holding arithmetic
 * types.
 * @tparam Pargs Types to be forwarded to recursive call of `save_varis`.
 * @param[in, out] dest Pointer to set of varis.
 * @param[in] x arithmetic or container of arithmetics.
 * @param[in] args Further arguments forwarded to recursive call.
 */
template <typename Arith, require_arithmetic_t<scalar_type_t<Arith>>* = nullptr,
          typename... Pargs>
inline vari** save_varis(vari** dest, Arith&& x, Pargs&&... args) {
  return save_varis(dest, std::forward<Pargs>(args)...);
}

/**
 * Ends the recursion by returning the input `vari**`
 * @param dest a pointer to pointers to varis.
 * @return the `dest` vari* pointer
 */
inline vari** save_varis(vari** dest) { return dest; }

}  // namespace math
}  // namespace stan

#endif
