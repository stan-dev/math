#ifndef STAN_MATH_PRIM_FUN_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/functor/function_gradients.hpp>

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
inline plain_type_t<EigMat> inverse(const EigMat& m) {
  ref_type_t<EigMat> m_ref = m;

  decltype(auto) val_fun = [](auto&& x) {
    check_square("inverse", "m", x);
    if (x.size() == 0) {
      return plain_type_t<decltype(x)>();
    }
    return x.inverse().eval();
  };
  auto rev_grad_fun = [](auto&& val, auto&& adj, auto&& shared_args, auto&& x) {
    return -multiply(multiply(transpose(val), adj), transpose(val));
  };
  auto fwd_grad_fun = [&](auto&& val, auto&& adj, auto&& shared_args, auto&& x) {
    return -multiply(multiply(val, adj), val);
  };
  auto shared_args_fun
      = [&](auto&& val, auto&& adj, auto&& x) {
        return std::tuple<double>(1.0);
      };
  return function_gradients(
      std::forward_as_tuple(m_ref),
      std::move(val_fun),
      std::forward<decltype(shared_args_fun)>(shared_args_fun),
      std::forward_as_tuple(rev_grad_fun),
      std::forward_as_tuple(fwd_grad_fun));
}

}  // namespace math
}  // namespace stan

#endif
