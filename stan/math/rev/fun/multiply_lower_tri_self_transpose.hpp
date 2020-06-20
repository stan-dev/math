#ifndef STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

namespace internal {

struct OpLLT : public AdjJacOp {
  int N_;
  int M_;
  double* l_mem_;

  template <std::size_t size, typename Derived>
  Eigen::MatrixXd operator()(const std::array<bool, size>& needs_adj,
			     const Eigen::MatrixBase<Derived>& L_arg) {
    const auto& L = L_arg.eval();
    
    N_ = L.rows();
    M_ = L.cols();

    l_mem_ = allocate_and_save(L.template triangularView<Eigen::Lower>());

    auto LT = map_matrix(l_mem_, N_, M_).transpose();

    return L.template triangularView<Eigen::Lower>() * LT;
  }

  template <std::size_t size, typename Derived>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::MatrixBase<Derived>& adj) {
    auto L = map_matrix(l_mem_, N_, M_);
 
    Eigen::MatrixXd adjL = (adj.transpose() + adj) * L;

    for(size_t j = 1; j < adjL.cols(); ++j)
      for(size_t i = 0; i < std::min(static_cast<size_t>(adjL.rows()), j); ++i)
	adjL(i, j) = 0.0;

    return std::make_tuple(adjL);
  }
};

}

template <typename T,
	  require_matrix2_t<T>* = nullptr,
	  require_any_var2_t<T>* = nullptr>
inline auto multiply_lower_tri_self_transpose(const T& L) {
  return adj_jac_apply<internal::OpLLT>(L);
}

}  // namespace math
}  // namespace stan
#endif
