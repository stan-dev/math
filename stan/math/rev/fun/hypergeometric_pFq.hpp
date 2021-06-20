#ifndef STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

namespace stan {
namespace math {

// var, var, var
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<is_var, Tp, Tq>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, const Tz& z) {
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

// var, double, var
template <typename Tp, typename Tq, typename Tz,
          require_eigen_vector_vt<is_var, Tp>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, Tq>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, const Tz& z) {
    arena_t<promote_scalar_t<var, Tp>> arena_p = p;
    arena_t<promote_scalar_t<double, Tq>> arena_q = q;
  return make_callback_var(hypergeometric_pFq(arena_p.val(), arena_q, z.val()),
                           [arena_p, arena_q, z](auto& vi) mutable {
                             vector_d grad_p(arena_p.size());
                             double grad_z;

                             grad_pFq_pz(grad_p, grad_z, arena_p.val(), arena_q, z.val());

                             arena_p.adj() += vi.adj() * grad_p;
                             z.adj() += vi.adj() * grad_z;
                           });
}

// var, double, double
template <typename Tp, typename Tq, typename Tz,
          require_eigen_vector_vt<is_var, Tp>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, Tq>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, const Tz& z) {
    arena_t<promote_scalar_t<var, Tp>> arena_p = p;
    arena_t<promote_scalar_t<double, Tq>> arena_q = q;
  return make_callback_var(hypergeometric_pFq(arena_p.val(), arena_q, z),
                           [arena_p, arena_q, z](auto& vi) mutable {
                             vector_d grad_p(arena_p.size());

                             grad_pFq_p(grad_p, arena_p.val(), arena_q, z);

                             arena_p.adj() += vi.adj() * grad_p;
                           });
}

// double, var, var
template <typename Tp, typename Tq, typename Tz,
          require_eigen_vector_vt<std::is_arithmetic, Tp>* = nullptr,
          require_eigen_vector_vt<is_var, Tq>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, const Tz& z) {
    arena_t<promote_scalar_t<double, Tp>> arena_p = p;
    arena_t<promote_scalar_t<var, Tq>> arena_q = q;
  return make_callback_var(hypergeometric_pFq(arena_p, arena_q.val(), z.val()),
                           [arena_p, arena_q, z](auto& vi) mutable {
                             vector_d grad_q(arena_q.size());
                             double grad_z;

                             grad_pFq_qz(grad_q, grad_z, arena_p, arena_q.val(), z.val());

                             arena_q.adj() += vi.adj() * grad_q;
                             z.adj() += vi.adj() * grad_z;
                           });
}

// double, var, double
template <typename Tp, typename Tq, typename Tz,
          require_eigen_vector_vt<std::is_arithmetic, Tp>* = nullptr,
          require_eigen_vector_vt<is_var, Tq>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, const Tz& z) {
    arena_t<promote_scalar_t<double, Tp>> arena_p = p;
    arena_t<promote_scalar_t<var, Tq>> arena_q = q;
  return make_callback_var(hypergeometric_pFq(arena_p, arena_q.val(), z),
                           [arena_p, arena_q, z](auto& vi) mutable {
                             vector_d grad_q(arena_q.size());

                             grad_pFq_q(grad_q, arena_p, arena_q.val(), z);

                             arena_q.adj() += vi.adj() * grad_q;
                           });
}

// double, double, var
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Tp, Tq>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, const Tz& z) {
    arena_t<promote_scalar_t<double, Tp>> arena_p = p;
    arena_t<promote_scalar_t<double, Tq>> arena_q = q;
  return make_callback_var(hypergeometric_pFq(arena_p, arena_q, z.val()),
                           [arena_p, arena_q, z](auto& vi) mutable {
                             double grad_z;

                             grad_pFq_z(grad_z, arena_p, arena_q, z.val());

                             z.adj() += vi.adj() * grad_z;
                           });
}

// var, var, double
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<is_var, Tp, Tq>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Tp& p, const Tq& q, const Tz& z) {
    arena_t<promote_scalar_t<var, Tp>> arena_p = p;
    arena_t<promote_scalar_t<var, Tq>> arena_q = q;
  return make_callback_var(hypergeometric_pFq(arena_p.val(), arena_q.val(), z),
                           [arena_p, arena_q, z](auto& vi) mutable {
                             vector_d grad_p(arena_p.size());
                             vector_d grad_q(arena_q.size());

                             grad_pFq_pq(grad_p, grad_q, arena_p.val(), arena_q.val(), z);

                             arena_p.adj() += vi.adj() * grad_p;
                             arena_q.adj() += vi.adj() * grad_q;
                           });
}

}  // namespace math
}  // namespace stan  
#endif
