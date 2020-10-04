#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr>
inline auto mdivide_left(const T1& A, const T2& B) {
  using ret_type = plain_type_t<decltype(A * B)>;
  
  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "B", B);

  if (A.size() == 0) {
    return ret_type(0, B.cols());
  }

  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);

  arena_t<promote_scalar_t<var, T1>> arena_A;
  arena_t<promote_scalar_t<var, T2>> arena_B;

  if (!is_constant<T1>::value) {
    arena_A = A_ref;
  }

  if (!is_constant<T2>::value) {
    arena_B = B_ref;
  }

  auto arena_A_val = to_arena(value_of(A_ref));
  arena_t<ret_type>
      res = arena_A_val.householderQr().solve(value_of(B_ref));

  reverse_pass_callback([arena_A, arena_B, arena_A_val, res]() mutable {
    promote_scalar_t<double, T2> adjB
      = arena_A_val.transpose().householderQr().solve(res.adj());

    if (!is_constant<T1>::value)
      arena_A.adj() += -adjB * res.val().transpose().eval();

    if (!is_constant<T2>::value)
      arena_B.adj() += adjB;
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
