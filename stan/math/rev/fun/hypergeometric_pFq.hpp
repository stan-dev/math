#ifndef STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_REV_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_v, vector_v, var)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<is_var, Ta, Tb>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<promote_scalar_t<var, Ta>> arena_a = a;
  arena_t<promote_scalar_t<var, Tb>> arena_b = b;
  return make_callback_var(
    hypergeometric_pFq(arena_a.val(), arena_b.val(), z.val()),
    [arena_a, arena_b, z](auto& vi) mutable {
      vector_d grad_a(arena_a.size());
      vector_d grad_b(arena_b.size());
      double grad_z;

      grad_pFq(grad_a, grad_b, grad_z, arena_a.val(), arena_b.val(), z.val());

      arena_a.adj() += vi.adj() * grad_a;
      arena_b.adj() += vi.adj() * grad_b;
      z.adj() += vi.adj() * grad_z;
    });
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_v, vector_d, var)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_eigen_vector_vt<is_var, Ta>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, Tb>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<promote_scalar_t<var, Ta>> arena_a = a;
  arena_t<promote_scalar_t<double, Tb>> arena_b = b;
  return make_callback_var(
    hypergeometric_pFq(arena_a.val(), arena_b, z.val()),
    [arena_a, arena_b, z](auto& vi) mutable {
      vector_d grad_a(arena_a.size());
      double grad_z;

      grad_pFq_pz(grad_a, grad_z, arena_a.val(), arena_b, z.val());

      arena_a.adj() += vi.adj() * grad_a;
      z.adj() += vi.adj() * grad_z;
    });
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_v, vector_d, double)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_eigen_vector_vt<is_var, Ta>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<promote_scalar_t<var, Ta>> arena_a = a;
  arena_t<promote_scalar_t<double, Tb>> arena_b = b;
  return make_callback_var(
    hypergeometric_pFq(arena_a.val(), arena_b, z),
    [arena_a, arena_b, z](auto& vi) mutable {
      vector_d grad_a(arena_a.size());

      grad_pFq_p(grad_a, arena_a.val(), arena_b, z);

      arena_a.adj() += vi.adj() * grad_a;
    });
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_v, var)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_eigen_vector_vt<std::is_arithmetic, Ta>* = nullptr,
          require_eigen_vector_vt<is_var, Tb>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<promote_scalar_t<double, Ta>> arena_a = a;
  arena_t<promote_scalar_t<var, Tb>> arena_b = b;
  return make_callback_var(
    hypergeometric_pFq(arena_a, arena_b.val(), z.val()),
    [arena_a, arena_b, z](auto& vi) mutable {
      vector_d grad_b(arena_b.size());
      double grad_z;

      grad_pFq_qz(grad_b, grad_z, arena_a, arena_b.val(), z.val());

      arena_b.adj() += vi.adj() * grad_b;
      z.adj() += vi.adj() * grad_z;
    });
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_v, double)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_eigen_vector_vt<std::is_arithmetic, Ta>* = nullptr,
          require_eigen_vector_vt<is_var, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<promote_scalar_t<double, Ta>> arena_a = a;
  arena_t<promote_scalar_t<var, Tb>> arena_b = b;
  return make_callback_var(
    hypergeometric_pFq(arena_a, arena_b.val(), z),
    [arena_a, arena_b, z](auto& vi) mutable {
      vector_d grad_b(arena_b.size());

      grad_pFq_q(grad_b, arena_a, arena_b.val(), z);

      arena_b.adj() += vi.adj() * grad_b;
    });
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_d, var)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Ta, Tb>* = nullptr,
          require_var_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<promote_scalar_t<double, Ta>> arena_a = a;
  arena_t<promote_scalar_t<double, Tb>> arena_b = b;
  return make_callback_var(hypergeometric_pFq(arena_a, arena_b, z.val()),
                          [arena_a, arena_b, z](auto& vi) mutable {
                            double grad_z;

                            grad_pFq_z(grad_z, arena_a, arena_b, z.val());

                            z.adj() += vi.adj() * grad_z;
                          });
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_v, vector_v, double)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<is_var, Ta, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline var hypergeometric_pFq(const Ta& a, const Tb& b, const Tz& z) {
  arena_t<promote_scalar_t<var, Ta>> arena_a = a;
  arena_t<promote_scalar_t<var, Tb>> arena_b = b;
  return make_callback_var(
    hypergeometric_pFq(arena_a.val(), arena_b.val(), z),
    [arena_a, arena_b, z](auto& vi) mutable {
      vector_d grad_a(arena_a.size());
      vector_d grad_b(arena_b.size());

      grad_pFq_pq(grad_a, grad_b, arena_a.val(), arena_b.val(), z);

      arena_a.adj() += vi.adj() * grad_a;
      arena_b.adj() += vi.adj() * grad_b;
    });
}

}  // namespace math
}  // namespace stan
#endif
