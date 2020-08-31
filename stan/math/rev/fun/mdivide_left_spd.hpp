#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_SPD_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T1, typename T2,
	  require_all_eigen_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr>
inline Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_spd(const T1& A, const T2& B) {
  check_symmetric("mdivide_left_spd", "A", A);
  check_multiplicable("mdivide_left_spd", "A", A, "B", B);

  using A_ref_t = ref_type_t<T1>;
  using B_ref_t = ref_type_t<T2>;

  A_ref_t A_ref = A;
  B_ref_t B_ref = B;
  
  check_not_nan("mdivide_left_spd", "A", A_ref);

  if (A.size() == 0) {
    return {0, B.cols()};
  }

  arena_matrix<promote_scalar_t<var, T1>> arena_A;
  arena_matrix<promote_scalar_t<var, T2>> arena_B;

  if (!is_constant<T1>::value) {
    arena_A = A_ref;
  }

  if (!is_constant<T2>::value) {
    arena_B = B_ref;
  }

  auto A_llt = value_of(A_ref).llt();
  check_pos_definite("mdivide_left_spd", "A", A_llt);

  arena_matrix<Eigen::MatrixXd> arena_A_llt(A.rows(), A.cols());
  arena_A_llt.template triangularView<Eigen::Lower>() = A_llt.matrixL();
  arena_A_llt.template triangularView<Eigen::Upper>() = arena_A_llt.transpose();  

  arena_matrix<Eigen::Matrix<var, T1::RowsAtCompileTime, T2::ColsAtCompileTime>>
    res = A_llt.solve(value_of(B_ref));

  reverse_pass_callback([arena_A, arena_B, arena_A_llt, res]() mutable {
    Eigen::Matrix<double, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
      adjB = res.adj();

    arena_A_llt.template triangularView<Eigen::Lower>().solveInPlace(adjB);
    arena_A_llt.template triangularView<Eigen::Upper>().solveInPlace(adjB);

    if(!is_constant<T1>::value)
      forward_as<promote_scalar_t<var, T1>>(arena_A).adj() +=
	-adjB * res.val().transpose().eval();

    if (!is_constant<T2>::value)
      forward_as<promote_scalar_t<var, T2>>(arena_B).adj() += adjB;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
