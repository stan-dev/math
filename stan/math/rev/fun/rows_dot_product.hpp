#ifndef STAN_MATH_REV_FUN_ROWS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_ROWS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {

struct RowsDotProductOp : public AdjJacOp {
  int N_;
  int M_;
  double* a_mem_;
  double* b_mem_;

  template <std::size_t size,
	    typename Derived1,
	    typename Derived2>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
			     const Eigen::MatrixBase<Derived1>& a_arg,
			     const Eigen::MatrixBase<Derived2>& b_arg) {
    const auto& a = a_arg.eval();
    const auto& b = b_arg.eval();
    
    N_ = a.rows();
    M_ = a.cols();

    if(needs_adj[0])
      b_mem_ = allocate_and_save(b);

    if(needs_adj[1])
      a_mem_ = allocate_and_save(a);

    Eigen::VectorXd out(N_);
    for(size_t n = 0; n < N_; ++n)
      out(n) = a.row(n).dot(b.row(n));
    
    return out;
  }

  template <std::size_t size, typename Derived>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::MatrixBase<Derived>& adj) {
    Eigen::MatrixXd adja;
    Eigen::MatrixXd adjb;

    if(needs_adj[0]) {
      auto b = map_matrix(b_mem_, N_, M_);
      adja.resize(N_, M_);

      for(size_t n = 0; n < N_; ++n)
	adja.row(n) = adj(n) * b.row(n);
    }

    if(needs_adj[1]) {
      auto a = map_matrix(a_mem_, N_, M_);
      adjb.resize(N_, M_);

      for(size_t n = 0; n < N_; ++n)
	adjb.row(n) = adj(n) * a.row(n);
    }

    return std::make_tuple(adja, adjb);
  }
};

}

template <typename Mat1, typename Mat2,
          require_matrix2_t<Mat1>* = nullptr,
          require_matrix2_t<Mat2>* = nullptr,
          require_any_var2_t<Mat1, Mat2>* = nullptr>
inline auto rows_dot_product(const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("rows_dot_product", "v1", v1, "v2", v2);

  return adj_jac_apply<internal::RowsDotProductOp>(v1, v2);
}

}  // namespace math
}  // namespace stan
#endif
