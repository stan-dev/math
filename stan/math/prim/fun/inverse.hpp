#ifndef STAN_MATH_PRIM_FUN_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/function_gradients_adj_jac.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not square.
 */
template <typename EigMat>
inline plain_type_t<EigMat>
inverse(const EigMat& m) {
  check_square("inverse", "m", m);
  if (m.size() == 0) {
    return {};
  }
  using par_t = partials_type_t<scalar_type_t<EigMat>>;
  decltype(auto) val_fun = [](auto&& x) { return inverse(x); };
  auto rev_grad_fun = [](auto&& val, auto&& adj, auto&& x) {
    decltype(auto) transpose_val = transpose(val);
    return -multiply(multiply(transpose_val, adj), transpose_val);
  };
  auto fwd_grad_fun = [](auto&& val, auto&& adj, auto&& x) {
    return -multiply(multiply(val, adj), val);
  };
  return function_gradients_adj_jac(std::forward_as_tuple(m),
                            std::move(val_fun),
                            std::forward_as_tuple(rev_grad_fun),
                            std::forward_as_tuple(fwd_grad_fun));
}

}  // namespace math
}  // namespace stan

#endif
