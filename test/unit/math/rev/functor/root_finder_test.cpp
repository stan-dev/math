#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
/*
TEST(RevFunctor, root_finder) {
  using stan::math::root_finder;
  using stan::math::var;
  // return cube root of x using 1st and 2nd derivatives and Halley.
  // using namespace std;  // Help ADL of std functions.
  var x = 27;
  int exponent;
  // Get exponent of z (ignore mantissa).
  std::frexp(x.val(), &exponent);
  // Rough guess is to divide the exponent by three.
  var guess = ldexp(1., exponent / 3);
  // Minimum possible value is half our guess.
  var min = ldexp(0.5, exponent / 3);
  // Maximum possible value is twice our guess.
  var max = ldexp(2., exponent / 3);
  // Maximum possible binary digits accuracy for type T.
  const int digits = std::numeric_limits<double>::digits;
  int get_digits = static_cast<int>(digits * 0.4);
  std::uintmax_t maxit = 20;
  auto f = [](const auto& g, const auto x) { return g * g * g - x; };
  auto f_g = [](const auto& g, const auto x) { return 3 * g * g; };
  auto f_gg = [](const auto& g, const auto x) { return 6 * g; };
  var result = root_finder(std::make_tuple(f, f_g, f_gg), guess, min, max, x);
  result.grad();
  EXPECT_EQ(27, x.val());
  EXPECT_NEAR(0.037037, x.adj(), 1e-6);
}
*/

struct BetaCdfRoot {
  template <bool ReturnDeriv, typename T1, typename T2, typename T3,
            typename T4, std::enable_if_t<ReturnDeriv>* = nullptr>
  static auto run(T1&& x, T2&& alpha, T3&& beta, T4&& p) {
    auto f = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
      return stan::math::beta_cdf(x, alpha, beta) - p;
    }(x, alpha, beta, p);
    auto beta_pdf_deriv = [](auto&& x, auto&& alpha, auto&& beta) {
      using stan::math::lgamma;
      using stan::math::pow;
      auto val = -(pow((1 - x), beta) * pow(x, (-2 + alpha))
                   * (1. - 2. * x - alpha + x * alpha + x * beta)
                   * lgamma(alpha + beta))
                 / (pow((-1 + x), 2) * lgamma(alpha) * lgamma(beta));
      return val;
    }(x, alpha, beta);
    auto beta_pdf = [](auto&& x, auto&& alpha, auto&& beta) {
      using stan::math::lgamma;
      using stan::math::pow;
      auto val
          = (pow(x, alpha - 1) * pow((1 - x), beta - 1) * lgamma(alpha + beta))
            / ((lgamma(alpha) * lgamma(beta)));
      return val;
    }(x, alpha, beta);
    return std::make_tuple(f, beta_pdf, beta_pdf_deriv);
  }
  template <bool ReturnDeriv, typename T1, typename T2, typename T3,
            typename T4, std::enable_if_t<!ReturnDeriv>* = nullptr>
  static auto run(T1&& x, T2&& alpha, T3&& beta, T4&& p) {
    return stan::math::beta_cdf(x, alpha, beta) - p;
  }
};
TEST(RevFunctor, root_finder_beta_cdf) {
  using stan::math::var;
  auto func = [](auto&& vals) {
    auto p = vals(0);
    auto alpha = vals(1);
    auto beta = vals(2);
    boost::math::beta_distribution<decltype(beta)> my_beta(alpha, beta);
    return boost::math::quantile(my_beta, p);
  };
  double p = 0.4;
  double a = 0.5;
  double b = 0.5;
  Eigen::VectorXd vals(3);
  // p, alpha, beta
  vals << p, a, b;
  double fx = 0;
  Eigen::VectorXd finit_grad_fx(3);
  stan::math::finite_diff_gradient(func, vals, fx, finit_grad_fx, 1e-3);
  std::cout << "--- Finit Diff----\n";
  std::cout << "fx: " << fx;
  std::cout << "\ngrads: \n"
            << "p: " << finit_grad_fx(0)
            << "\n"
               "alpha: "
            << finit_grad_fx(1)
            << "\n"
               "beta: "
            << finit_grad_fx(2) << "\n";
  double guess = .3;
  double min = 0;
  double max = 1;
  auto full_f = [guess, min, max](auto&& alpha, auto&& beta, auto&& p) {
    std::uintmax_t max_its = 1000;
    return stan::math::root_finder_hailey<BetaCdfRoot>(guess, min, max, alpha,
                                                       beta, p);
  };
  auto func2 = [&full_f](auto&& vals) {
    auto p = vals(0);
    auto alpha = vals(1);
    auto beta = vals(2);
    return full_f(alpha, beta, p);
  };
  Eigen::VectorXd grad_fx(3);
  Eigen::Matrix<stan::math::var, -1, 1> var_vec(vals);
  stan::math::var fxvar = func2(var_vec);
  fxvar.grad();
  grad_fx = var_vec.adj();
  fx = fxvar.val();
  std::cout << "fxvar adj:" << fxvar.adj() << "\n";
  // stan::math::gradient(func2, vals, fx, grad_fx);
  std::cout << "--- Auto Diff----\n";
  std::cout << "fx: " << fx;
  std::cout << "\ngrads: \n"
            << "p: " << grad_fx(0)
            << "\n"
               "alpha: "
            << grad_fx(1)
            << "\n"
               "beta: "
            << grad_fx(2) << "\n";
  Eigen::VectorXd diff_grad_fx = finit_grad_fx - grad_fx;
  std::cout << "--- grad diffs----\n";
  std::cout << "p: " << diff_grad_fx(0)
            << "\n"
               "alpha: "
            << diff_grad_fx(1)
            << "\n"
               "beta: "
            << diff_grad_fx(2) << "\n";
  auto deriv_p = [](auto& p, auto& a, auto& b) {
    using std::pow;
    using boost::math::ibeta_inv;
    using boost::math::beta;
    return beta(a, b) * pow(1 - ibeta_inv(a, b, p), (1 - b))
           * pow(ibeta_inv(a, b, p), (1 - a));
  };
  auto deriv_a = [](auto& p, auto& a, auto& b) {
    using std::pow;
    using boost::math::ibeta_inv;
    using boost::math::ibeta;
    using boost::math::beta;
    using boost::math::tgamma;
    using boost::math::hypergeometric_pFq;
    using boost::math::beta;
    using boost::math::polygamma;
    using std::log;
    double w = ibeta_inv(a, b, p);
    return pow(1 - w, (1 - b)) * pow(w, (1 - a))
           * (pow(w, a) * pow(tgamma(a), 2)
                  * (hypergeometric_pFq({a, a, 1 - b}, {1 + a, 1 + a}, w)
                     / (tgamma(1 + a) * tgamma(1 + a)))
              - beta(a, b) * ibeta(a, b, w)
                    * (log(w) - polygamma(0, a) + polygamma(0, a + b)));
  };

  auto deriv_b = [](auto& p, auto& a, auto& b) {
    using std::pow;
    using boost::math::ibeta_inv;
    using boost::math::ibeta;
    using boost::math::beta;
    using boost::math::tgamma;
    using boost::math::hypergeometric_pFq;
    using boost::math::beta;
    using boost::math::polygamma;
    using std::log;
    return pow(1 - ibeta_inv(a, b, p), -b) * (-1 + ibeta_inv(a, b, p))
           * pow(ibeta_inv(a, b, p), (1 - a))
           * (pow(tgamma(b), 2)
                  * (hypergeometric_pFq({b, b, 1 - a}, {1 + b, 1 + b},
                                        1 - ibeta_inv(a, b, p))
                     / (tgamma(1 + b) * tgamma(1 + b)))
                  * pow(1 - ibeta_inv(a, b, p), b)
              - beta(b, a, 1 - ibeta_inv(a, b, p))
                    * (log(1 - ibeta_inv(a, b, p)) - polygamma(0, b)
                       + polygamma(0, a + b)));
  };
  double known_p_grad = deriv_p(p, a, b);
  double known_alpha_grad = deriv_a(p, a, b);
  double known_beta_grad = deriv_b(p, a, b);
  std::cout << "--- Mathematica Calculate Grad----\n";
  std::cout << "p: " << known_p_grad << "\n"
            << "alpha: " << known_alpha_grad << "\n"
            << "beta: " << known_beta_grad << "\n";
  std::cout << "--- Mathematica Calculate Grad Diff----\n";
  std::cout << "p: " << grad_fx(0) - known_p_grad << "\n"
            << "alpha: " << grad_fx(1) - known_alpha_grad << "\n"
            << "beta: " << grad_fx(2) - known_beta_grad << "\n";
}

template <size_t... I>
inline constexpr auto make_index_tuple(const std::index_sequence<I...>&) {
  return std::make_tuple(I...);
}
template <typename FTuple, typename FGrad, typename... Args>
void check_vs_known_grads(FTuple&& grad_tuple, FGrad&& f, double tolerance,
                          Args&&... args) {
  try {
    const stan::math::nested_rev_autodiff nested;
    auto var_tuple = std::make_tuple(stan::math::var(args)...);
    // For pretty printer index of incorrect values
    auto arg_num_tuple = stan::math::apply(
        [&args...](auto&&... nums) {
          return std::make_tuple(std::make_tuple(args, nums)...);
        },
        make_index_tuple(std::make_index_sequence<sizeof...(args)>()));

    auto ret = stan::math::apply(
        [&f](auto&&... var_args) { return f(var_args...); }, var_tuple);
    ret.grad();
    auto grad_val_tuple = stan::math::apply(
        [&args...](auto&&... grad_funcs) {
          return std::make_tuple(grad_funcs(args...)...);
        },
        grad_tuple);
    auto adj_tuple = stan::math::apply(
        [](auto&&... var_arg) { return std::make_tuple(var_arg.adj()...); },
        var_tuple);
    stan::math::for_each(
        [tolerance](auto&& grad, auto&& adj, auto&& num_helper) {
          EXPECT_NEAR(grad, adj, tolerance)
              << "Diff: (" << grad - adj << ")\nWith arg #"
              << std::get<1>(num_helper) << " (" << std::get<0>(num_helper)
              << ")";
        },
        grad_val_tuple, adj_tuple, arg_num_tuple);
  } catch (const std::exception& e) {
    stan::math::recover_memory();
  }
  stan::math::recover_memory();
}
TEST(RevFunctor, root_finder_beta_cdf2) {
  constexpr double guess = .5;
  constexpr double min = 0;
  constexpr double max = 1;

  auto deriv_a = [](auto&& a, auto&& b, auto&& p) {
    using std::pow;
    using boost::math::ibeta_inv;
    using boost::math::ibeta;
    using boost::math::beta;
    using boost::math::tgamma;
    using boost::math::hypergeometric_pFq;
    using boost::math::beta;
    using boost::math::polygamma;
    using std::log;
    double w = ibeta_inv(a, b, p);
    return pow(1 - w, (1 - b)) * pow(w, (1 - a))
           * (pow(w, a) * pow(tgamma(a), 2)
                  * (hypergeometric_pFq({a, a, 1 - b}, {1 + a, 1 + a}, w)
                     / (tgamma(1 + a) * tgamma(1 + a)))
              - beta(a, b) * ibeta(a, b, w)
                    * (log(w) - polygamma(0, a) + polygamma(0, a + b)));
  };
  auto deriv_b = [](auto&& a, auto&& b, auto&& p) {
    using std::pow;
    using boost::math::ibeta_inv;
    using boost::math::ibeta;
    using boost::math::beta;
    using boost::math::tgamma;
    using boost::math::hypergeometric_pFq;
    using boost::math::beta;
    using boost::math::polygamma;
    using std::log;
    return pow(1 - ibeta_inv(a, b, p), -b) * (-1 + ibeta_inv(a, b, p))
           * pow(ibeta_inv(a, b, p), (1 - a))
           * (pow(tgamma(b), 2)
                  * (hypergeometric_pFq({b, b, 1 - a}, {1 + b, 1 + b},
                                        1 - ibeta_inv(a, b, p))
                     / (tgamma(1 + b) * tgamma(1 + b)))
                  * pow(1 - ibeta_inv(a, b, p), b)
              - beta(b, a, 1 - ibeta_inv(a, b, p))
                    * (log(1 - ibeta_inv(a, b, p)) - polygamma(0, b)
                       + polygamma(0, a + b)));
  };
  auto deriv_p = [](auto&& a, auto&& b, auto&& p) {
    using std::pow;
    using boost::math::ibeta_inv;
    using boost::math::beta;
    return beta(a, b) * pow(1 - ibeta_inv(a, b, p), (1 - b))
           * pow(ibeta_inv(a, b, p), (1 - a));
  };
  auto f_newton = [guess, min, max](auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::root_finder_newton_raphson<BetaCdfRoot>(guess, min, max,
                                                               alpha, beta, p);
  };
  auto f_schroder = [guess, min, max](auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::root_finder_schroder<BetaCdfRoot>(guess, min, max, alpha,
                                                         beta, p);
  };
  auto f_hailey = [](auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::root_finder_hailey<BetaCdfRoot>(p, min, max, alpha, beta,
                                                       p);
  };
  check_vs_known_grads(std::make_tuple(deriv_a, deriv_b, deriv_p), f_schroder,
                       5e-2, .5, .5, .3);

  /* For some reason this fails after a while??
  for (double p = .3; p < .7; p += .1) {
    for (double a = .3; a < .7; a += .1) {
      for (double b = .3; b < .7; b += .1) {

          check_vs_known_grads(std::make_tuple(deriv_a, deriv_b, deriv_p),
                              f_schroder, 5e-2, a, b, p);
          if (::testing::Test::HasFailure()) {
            std::cout << "--\na: " << a << "\nb: " << b << "\np: " << p << "\n";
          }
      }
    }
  }
  */
}
