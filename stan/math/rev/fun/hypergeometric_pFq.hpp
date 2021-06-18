#ifndef STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

namespace stan {
namespace math {

template <typename Tp, typename Tq,
          require_all_eigen_vector_st<is_var, Tp, Tq>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, var z) {
    arena_t<promote_scalar_t<var, Tp>> arena_p = p;
    arena_t<promote_scalar_t<var, Tq>> arena_q = q;
  return make_callback_var(hypergeometric_pFq(arena_p.val(), arena_q.val(), z.val()),
                           [arena_p, arena_q, z](auto& vi) mutable {
                             vector_d grad_p(arena_p.size());
                             vector_d grad_q(arena_q.size());
                             double grad_z;

                             grad_pFq(grad_p, grad_q, grad_z, arena_p.val(), arena_q.val(), z.val());

                             arena_p.adj() += vi.adj() * grad_p;
                             arena_q.adj() += vi.adj() * grad_q;
                             z.adj() += vi.adj() * grad_z;
                           });
}

}  // namespace math
}  // namespace stan  
#endif
