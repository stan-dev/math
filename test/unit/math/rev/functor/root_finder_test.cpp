#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions/beta.hpp>
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
TEST(RevFunctor, root_finder_beta_cdf) {
  using stan::math::var;
  auto func = [](auto&& vals) {
    auto p = vals(0);
    auto alpha = vals(1);
    auto beta = vals(2);
    boost::math::beta_distribution<decltype(beta)> my_beta(alpha, beta);
    return boost::math::quantile(my_beta, p);
  };
  Eigen::VectorXd vals(3);
  // p, alpha, beta
  vals << 0.4, .5, .5;
  double fx = 0;
  Eigen::VectorXd finit_grad_fx(3);
  stan::math::finite_diff_gradient(func, vals, fx, finit_grad_fx, 1e-14);
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

  auto beta_pdf_deriv = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    using stan::math::lgamma;
    using stan::math::pow;
    auto val = -(pow((1 - x), beta) * pow(x, (-2 + alpha))
                 * (1. - 2. * x - alpha + x * alpha + x * beta)
                 * lgamma(alpha + beta))
               / (pow((-1 + x), 2) * lgamma(alpha) * lgamma(beta));
    return val;
  };
  auto beta_pdf = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    using stan::math::lgamma;
    using stan::math::pow;
    auto val
        = (pow(x, alpha - 1) * pow((1 - x), beta - 1) * lgamma(alpha + beta))
          / ((lgamma(alpha) * lgamma(beta)));
    return val;
  };

  auto f = [](auto&& x, auto&& alpha, auto&& beta, auto&& p) {
    return stan::math::beta_cdf(x, alpha, beta) - p;
  };
  double guess = .3;
  double min = 0;
  double max = 1;
  auto full_f = [&f, &beta_pdf_deriv, &beta_pdf, guess, min, max](
                    auto&& alpha, auto&& beta, auto&& p) {
    std::uintmax_t max_its = 1000;
    return stan::math::root_finder_hailey(
        std::make_tuple(f, beta_pdf, beta_pdf_deriv), guess, min, max,
        alpha, beta, p);
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
  //stan::math::gradient(func2, vals, fx, grad_fx);
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

  double known_p_grad = 1.493914820513846;
  double known_alpha_grad = 0.948540678312004;
  double known_beta_grad = -0.7464041618483364;
  std::cout << "--- Auto Diff----\n";
  std::cout << "p: " << grad_fx(0) - known_p_grad << "\n" <<
   "alpha: " << grad_fx(1) - known_alpha_grad << "\n" <<
   "beta: " << grad_fx(2) - known_beta_grad << "\n";
}
