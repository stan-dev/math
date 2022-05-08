#ifndef STAN_MATH_PRIM_FUN_INV_PHI_LOG_HPP
#define STAN_MATH_PRIM_FUN_INV_PHI_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
/**
 * The inverse of the unit normal cumulative distribution function.
 *
 * @param p argument between 0 and 1 inclusive
 * @return Real value of the inverse cdf for the standard normal distribution.
 */
inline long double inv_Phi_log_lambda(long double log_p) {
  check_less_or_equal("inv_Phi_log", "Probability variable", log_p, 0.);

  if (log_p == NEGATIVE_INFTY) {
    return NEGATIVE_INFTY;
  }
  if (log_p == 0) {
    return INFTY;
  }

  static const long double log_a[8]
      = {1.2199838032983212361078285710532819636747869949455269702015,
         4.8914137334471356660028877181679071443394481346972814858613,
         7.5865960847956080883082295570387898324568001714816000563527,
         9.5274618535358388901978394767678575968439856675022781953217,
         10.734698580862359193568558400271381214615993152487675386540,
         11.116406781896242115649156409480902760271809094853583122134,
         10.417226196842595017693249813096834274194835611478386401665,
         7.8276718012189362623726348876474997383871412260093333891226};
  static const long double log_b[8]
      = {0.,
         3.7451021830139207383037397051332758147125149878969524050480,
         6.5326064640478618289295562233751573956126769356556449255276,
         8.5930788436817044087736811921133610029011365371840851909132,
         9.9624069236663077425466579299794622162198838778518543394569,
         10.579180688621286666206645046714446503611903777605610455541,
         10.265665328832871792104231874798798249892953877879038702305,
         8.5614962136628454945793571514493658608756478643490687605256};
  static const long double log_c[8]
      = {0.3530744474482423474253411971863875773110416729361317536072,
         1.5326298343683388653365882077386965677836566047627135895898,
         1.7525849400614634832553515573234472525649040492075733723068,
         1.2941374937060454491309878479045845831522022489649421832255,
         0.2393776640901312565674683205114411320646903467410118642625,
         -1.419724057885092049733742418443070998814921203853709341365,
         -3.784340465764968143046047427817485545086795726249017944579,
         -7.163234779359426533991313188096094249706979717037172942753};
  static const long double log_d[8]
      = {0.,
         0.7193954734947205047245123628276879846369669216368876768165,
         0.5166395879845317615211151540976544006541209051453173240568,
         -0.371400933927844351706587139965567469873214555879242076139,
         -1.909840708457214045519577504534589395147865025037188124762,
         -4.186547581055928432171715090347894581668539504453380854814,
         -7.509976771225415296023136809870347991212091151984611110563,
         -20.67376157385924942206085235134755772380694892116629078845};
  static const long double log_e[8]
      = {1.8958048169567149828368030301222133925312189249850017925505,
         1.6981417567726154500599548596591401682246697139577459172901,
         0.5793212339927351466920742848614918786577535875989158702308,
         -1.215503791936417740743353960363053177354231954275848136327,
         -3.629396584023968393125053306667390321132075537098856528171,
         -6.690500273261249798544839142323179876784895297437975706380,
         -10.51540298415323795960908434918974531624763999447252923447,
         -15.41979457491781577083604431696521397005610506703480596500};
  static const long double log_f[8]
      = {0.,
         -0.511105318617135869400864118070037324958858441968973309655,
         -1.988286302259815913266915581038346510957576044950301014734,
         -4.208049039384857940412998843227725455979622001916832277710,
         -7.147448611626374901826054326105044483668972676628832825814,
         -10.89973190740069813377286491971487048164963921352208678904,
         -15.76637472711685566869394667681439515037736894877363436979,
         -33.82373901099482547859307148234549286520305771584510398829};

  long double val;
  long double log_q = log_p <= LOG_HALF
                          ? log_diff_exp(0, log_sum_exp(log_p, LOG_HALF))
                          : log_diff_exp(log_p, LOG_HALF);
  int log_q_sign = log_p <= LOG_HALF ? -1 : 1;

  if (log_q <= -0.855666110057720222602921079727590808163898744763538696977) {
    long double log_r = log_diff_exp(
        -1.711332220115440445205842159455181616327797489527077393954,
        2 * log_q);
    long double log_agg_a = log_sum_exp(log_a[7] + log_r, log_a[6]);
    long double log_agg_b = log_sum_exp(log_b[7] + log_r, log_b[6]);

    for (int i = 0; i < 6; i++) {
      log_agg_a = log_sum_exp(log_agg_a + log_r, log_a[5 - i]);
      log_agg_b = log_sum_exp(log_agg_b + log_r, log_b[5 - i]);
    }

    return log_q_sign * exp(log_q + log_agg_a - log_agg_b);
  } else {
    long double log_r = log_q_sign == -1 ? log_p : log_diff_exp(0., log_p);

    if (stan::math::is_inf(log_r)) {
      return 0;
    }

    log_r = log(sqrt(-log_r));

    if (log_r
        <= 1.6094379124341003746007593332261876395256013542685177219126478914) {
      log_r = log_diff_exp(
          log_r, 0.4700036292457355536509370311483420647008990488122480404493);
      long double log_agg_c = log_sum_exp(log_c[7] + log_r, log_c[6]);
      long double log_agg_d = log_sum_exp(log_d[7] + log_r, log_d[6]);

      for (int i = 0; i < 6; i++) {
        log_agg_c = log_sum_exp(log_agg_c + log_r, log_c[5 - i]);
        log_agg_d = log_sum_exp(log_agg_d + log_r, log_d[5 - i]);
      }
      val = exp(log_agg_c - log_agg_d);
    } else {
      log_r = log_diff_exp(
          log_r,
          1.6094379124341003746007593332261876395256013542685177219126478914);
      long double log_agg_e = log_sum_exp(log_e[7] + log_r, log_e[6]);
      long double log_agg_f = log_sum_exp(log_f[7] + log_r, log_f[6]);

      for (int i = 0; i < 6; i++) {
        log_agg_e = log_sum_exp(log_agg_e + log_r, log_e[5 - i]);
        log_agg_f = log_sum_exp(log_agg_f + log_r, log_f[5 - i]);
      }
      val = exp(log_agg_e - log_agg_f);
    }
    if (log_q_sign == -1)
      return -val;
  }
  return val;
}
}  // namespace internal

/**
 * Return the value of the inverse standard normal cumulative distribution
 * function at the specified argument.
 *
 * The precision is at or better than 1.5e-15 for values between 0.0000001 he
 * largest integer that protects against floating point errors for the inv_Phi
 * function. The value was found by finding the largest integer that passed the
 * unit tests for accuracy when the input into inv_Phi is near 1.
 *
 * @param p argument between 0 and 1 inclusive
 * @return real value of the inverse cdf for the standard normal distribution
 */
inline long double inv_Phi_log(long double log_p) {
  return internal::inv_Phi_log_lambda(log_p);
}

/**
 * Structure to wrap inv_Phi_log() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable in range [0, 1]
 * @return Inverse unit normal CDF of x.
 * @throw std::domain_error if x is not between 0 and 1.
 */
struct inv_Phi_log_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_Phi_log(x);
  }
};

/**
 * Vectorized version of inv_Phi_log().
 *
 * @tparam T type of container
 * @param x variables in range [0, 1]
 * @return Inverse unit normal CDF of each value in x.
 * @throw std::domain_error if any value is not between 0 and 1.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto inv_Phi_log(const T& x) {
  return apply_scalar_unary<inv_Phi_log_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
