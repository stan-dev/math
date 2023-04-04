#ifndef STAN_MATH_PRIM_PROB_WIENER_FULL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER_FULL_LPDF_HPP

#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob/wiener5_lpdf.hpp>
#include <stan/math/prim/prob/wiener7_lpdf.hpp>

namespace stan {
namespace math {
namespace internal {

template <bool propto, typename T_y, typename T_a, typename T_t0, typename T_w,
          typename T_v, typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>
wiener_full_prec_impl_lpdf(const char* function_name, const T_y& y,
                           const T_a& a, const T_t0& t0, const T_w& w,
                           const T_v& v, const T_sv& sv, const T_sw& sw,
                           const T_st0& st0, const double& prec) {
  using T_sw_ref = ref_type_t<T_sw>;
  using T_st0_ref = ref_type_t<T_st0>;

  check_consistent_size(function_name,
                        "Inter-trial variability in A-priori bias", sw, 1);
  check_consistent_size(function_name,
                        "Inter-trial variability in Nondecision time", st0, 1);

  T_sw_ref sw_ref = sw;
  T_st0_ref st0_ref = st0;

  size_t N = max_size(y, a, t0, w, v, sv, sw, st0);
  if (!N) {
    return 0;
  }

  scalar_seq_view<T_sw_ref> sw_vec(sw_ref);
  scalar_seq_view<T_st0_ref> st0_vec(st0_ref);

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    // Calculate 4-parameter model without inter-trial variabilities (if
    // sv_vec[i] == 0) or 5-parameter model with inter-trial variability in
    // drift rate (if sv_vec[i] != 0)
    if (sw_vec[i] == 0 && st0_vec[i] == 0) {
      return wiener5_lpdf<propto>(y, a, t0, w, v, sv, prec);
      // Calculate 6-, or 7-parameter model
    } else {
      return wiener7_lpdf<propto>(y, a, t0, w, v, sv, sw, st0, prec);
    }
  }
  return 0;
}
//-----------------------------------------------

}  // namespace internal

/** \ingroup prob_dists
 * The log of the first passage time density function for a (Wiener)
 * drift diffusion model with up to 7 parameters, where
 * \f$y\in \mathbb{R}_{+}\f$ is the reacion time, \f$a \in \mathbb{R}_{+}\f$
 * the boundary separation, \f$t_0 \in \mathbb{R}_{\geq 0}\f$ the non-decision
 time,
 * \f$w \in (0, 1)\f$ the relative starting point (aka a-priori bias),
 * \f$v \in \mathbb{R}\f$ the drifte rate, \f$s_v \in
 * \mathbb{R}_{\geq 0}\f$ the inter-trial variability of the drift rate,
 * \f$s_w \in [0, 1)\f$ the inter-trial  variability of the relative starting
 * point, and  \f$s_{t_0} \in \mathbb{R}_{\geq 0}\f$ the inter-trial variability
 * of
 * the non-decision time.
 *
 *
 \f{eqnarray*}{
 y &\sim& \text{wiener_full}(a,t_0,w,v,s_v,s_w,s_{t_0}) \\
 \log(p(y|a,v,w,t_0,s_v,s_w,s_{t_0})) &=& \log(\frac{1}{s_{t_0}}
 \int_{t_0}^{t_o + s_{t_0}} \frac{1}{s_{w}}\int_{w -0.5s_w}^{w + 0.5s_{w}}
 \int_{-\infty}^{\infty} p_3(y-\tau_0|a,\nu,\omega) \\
 &&\times\frac{1}{\sqrt{2\pi s_\nu^2}}
 \mathbb{e}^{-\frac{(\nu-v)^2}{2s_\nu^2}} \ d\nu \ d\omega \ d\tau_0) \\
 &=& \log(\frac{1}{s_{t_0}}
 \int_{t_0}^{t_o + s_{t_0}} \frac{1}{s_{w}}\int_{w -0.5s_w}^{w + 0.5s_{w}}
 M \times p_3(y-\tau_0|a,v,\omega) \ d\omega \ d\tau_0),
 \f}
 * where \f$M\f$ and \f$p_3()\f$ are defined, by using \f$t:=y-\tau_0\f$, as
 \f{eqnarray*}{
 M &:=& \frac{1}{\sqrt{1+s^2_v t}}
 \mathbb{e}^{av\omega+\frac{v^2t}{2}+\frac{s^2_v a^2
 \omega^2-2av\omega-v^2t}{2(1+s^2_vt)}} \text{ and} \\ p_3(t|a,v,w) &:=&
 \frac{1}{a^2} \mathbb{e}^{-a v w -\frac{v^2t}{2}} f(\frac{t}{a^2}|0,1,w), \f}
 * where  \f$f(t^*=\frac{t}{a^2}|0,1,w)\f$ has two forms
 \f{eqnarray*}{
 f_l(t^*|0,1,w) &=& \sum_{k=1}^{\infty} k\pi \mathbb{e}^{-\frac{k^2\pi^2t^*}{2}}
 \sin{(k \pi w)}\text{ and} \\
 f_s(t^*|0,1,w) &=& \sum_{k=-\infty}^{\infty} \frac{1}{\sqrt{2\pi (t^*)^3}}
 (w+2k) \mathbb{e}^{-\frac{(w+2k)^2}{2t^*}}, \f}
 * which are selected depending on the number of components \f$k\f$, needed to
 * guarantee a certain precision.
 *
 * Note that the parametrization for non-decision time and relative starting
 point is as follows:
 * \c t0 is the lower bound of the variability interval;
 * \c w is the mean of the variability interval.
 *
 * See \b Details below for more details on how to use \c wiener_full_lpdf().
 *
 * @tparam T_y type of scalar
 * @tparam T_a type of boundary
 * @tparam T_t0 type of non-decision time
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 * @tparam T_sv type of inter-trial variability of drift rate
 * @tparam T_sw type of inter-trial variability of relative starting point
 * @tparam T_st0 type of inter-trial variability of non-decision time
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param t0 The non-decision time
 * @param w The relative starting point
 * @param v The drift rate
 * @param sv The inter-trial variability of the drift rate
 * @param sw The inter-trial variability of the relative starting point
 * @param st0 The inter-trial variability of the non-decision time
 * @return The log of the Wiener first passage time density with
 *  the specified arguments for upper boundary responses
 * @throw std::domain_error if non-decision time \c t0 is greater than reaction
 time \c y.
 * @throw std::domain_error if \c 1-sw/2 is smaller than or equal to \c w.
 * @throw std::domain_error if \c sw/2 is larger than or equal to \c w.
 *
 *
 *
 * **Details**
 *
 * The function can be called by
 @code
 target += wiener_full_lpdf(y, a, t0, w, v, sv, sw, st0);
 @endcode
 * or
 @code
 y ~ wiener_full(a, t0, w, v, sv, sw, st0);
 @endcode
 *
 * By default \c wiener_full_lpdf() gives the log of the
 * Wiener first-passage time probability density function for the \e upper
 response
 * boundary. To use the \e lower response boundary \c v and \c w must be changed
 to
 * \c -v and \c (1-w), respectively.
 *
 * \c sv, \c sw, \c st0, and \c t0 can be set to zero, indicating no inter-trial
 * variability in \f$v\f$, no inter-trial variability in \f$w\f$, no
 * inter-trial variability in \f$t_0\f$, and no non-decision time, respectively.
 * If \c t0 is zero, \c st0 must be zero as well. For example, when no
 inter-trial
 * variability for the relative starting point is needed one can write something
 like:
 @code
 target += wiener_full_lpdf(y, a, t0, w, v, sv, 0, st0)
 @endcode
 * If no inter-trial variability is needed at all one can write something like:
 @code
 target += wiener_full_lpdf(y, a, t0, w, v, 0, 0, 0)
 @endcode
 * If for some reason no non-decision time is assumed one can write something
 like:
 @code
 target += wiener_full_lpdf(y, a, 0, w, v, sv, sw, 0)
 @endcode
 *
 * To also control the precision in the estimation of the partial derivatives,
 * call the function \c wiener_full_prec_lpdf(), analogously:
 @code
 target += wiener_full__prec_lpdf(y, a, t0, w, v, sv, sw, st0, precision);
 @endcode
 *
 *
 * **References**
 * - Blurton, S. P., Kesselmeier, M., & Gondan, M. (2017). The first-passage
 * time distribution for the diffusion model with variable drift.
 * *Journal of Mathematical Psychology, 76*, 7–12.
 * https://doi.org/10.1016/j.jmp.2016.11.003
 * - Foster, K., & Singmann, H. (2021). Another Approximation of the
 First-Passage
 * Time Densities for the Ratcliff Diffusion Decision Model.
 * *arXiv preprint arXiv:2104.01902*
 * - Gondan, M., Blurton, S. P., & Kesselmeier, M. (2014). Even faster and even
 * more accurate first-passage time densities and distributions for the Wiener
 * diffusion model. *Journal of Mathematical Psychology, 60*, 20–22.
 * https://doi.org/10.1016/j.jmp.2014.05.002
 * - Hartmann, R., & Klauer, K. C. (2021). Partial derivatives for the
 * first-passage time distribution in Wiener diffusion models.
 * *Journal of Mathematical Psychology, 103*, 102550.
 * https://doi.org/10.1016/j.jmp.2021.102550
 * - Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate calculations for
 * first-passage times in Wiener diffusion models.
 * *Journal of Mathematical Psychology, 53*(4), 222–230.
 * https://doi.org/10.1016/j.jmp.2009.02.003
 */

template <bool propto, typename T_y, typename T_a, typename T_t0, typename T_w,
          typename T_v, typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>
wiener_full_lpdf(const T_y& y, const T_a& a, const T_t0& t0, const T_w& w,
                 const T_v& v, const T_sv& sv, const T_sw& sw,
                 const T_st0& st0) {
  double precision = 1e-4;

  return internal::wiener_full_prec_impl_lpdf<propto, T_y, T_a, T_t0, T_w, T_v,
                                              T_sv, T_sw, T_st0>(
      "wiener_full_lpdf", y, a, t0, w, v, sv, sw, st0, precision);
}

template <typename T_y, typename T_a, typename T_t0, typename T_w, typename T_v,
          typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>
wiener_full_lpdf(const T_y& y, const T_a& a, const T_t0& t0, const T_w& w,
                 const T_v& v, const T_sv& sv, const T_sw& sw,
                 const T_st0& st0) {
  double precision = 1e-4;

  return internal::wiener_full_prec_impl_lpdf<false>(
      "wiener_full_lpdf", y, a, t0, w, v, sv, sw, st0, precision);
}

/** \ingroup prob_dists
 * The log of the first passage time density function for a (Wiener)
 * drift diffusion model with up to 7 parameters with the option
 * to control for the precision in the estimation of the partial derivatives.
 *
 * For \b Details see \c wiener_full_lpdf(). The usage and behavior of these
 * functions are the same excpet of the control over the precision.
 */

template <bool propto, typename T_y, typename T_a, typename T_t0, typename T_w,
          typename T_v, typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>
wiener_full_prec_lpdf(const T_y& y, const T_a& a, const T_t0& t0, const T_w& w,
                      const T_v& v, const T_sv& sv, const T_sw& sw,
                      const T_st0& st0, const double& prec) {
  return internal::wiener_full_prec_impl_lpdf<propto, T_y, T_a, T_t0, T_w, T_v,
                                              T_sv, T_sw, T_st0>(
      "wiener_full_prec_lpdf", y, a, t0, w, v, sv, sw, st0, prec);
}

template <typename T_y, typename T_a, typename T_t0, typename T_w, typename T_v,
          typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>
wiener_full_prec_lpdf(const T_y& y, const T_a& a, const T_t0& t0, const T_w& w,
                      const T_v& v, const T_sv& sv, const T_sw& sw,
                      const T_st0& st0, const double& prec) {
  return internal::wiener_full_prec_impl_lpdf<false>(
      "wiener_full_prec_lpdf", y, a, t0, w, v, sv, sw, st0, prec);
}

}  // namespace math
}  // namespace stan
#endif
