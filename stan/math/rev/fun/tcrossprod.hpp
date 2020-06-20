#ifndef STAN_MATH_REV_FUN_TCROSSPROD_HPP
#define STAN_MATH_REV_FUN_TCROSSPROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/Eigen_NumTraits.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

namespace internal {

struct TCrossProdOp : public AdjJacOp {
  int N_;
  int M_;
  double* x_mem_;

  template <std::size_t size, typename Derived>
  Eigen::MatrixXd operator()(const std::array<bool, size>& needs_adj,
			     const Eigen::MatrixBase<Derived>& x_arg) {
    const auto& x = x_arg.eval();
    
    N_ = x.rows();
    M_ = x.cols();

    x_mem_ = allocate_and_save(x);

    return x * x.transpose();
  }

  template <std::size_t size, typename Derived>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::MatrixBase<Derived>& adj) {
    auto x = map_matrix(x_mem_, N_, M_);
 
    Eigen::MatrixXd adjx = (adj.transpose() + adj) * x;

    return std::make_tuple(adjx);
  }
};

}
/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T Type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T,
	  require_matrix2_t<T>* = nullptr,
	  require_any_var2_t<T>* = nullptr>
inline auto tcrossprod(const T& M) {
  return adj_jac_apply<internal::TCrossProdOp>(M);
}

}  // namespace math
}  // namespace stan
#endif
