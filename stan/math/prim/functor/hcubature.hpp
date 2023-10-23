// For the code in hcubature.hpp, the package JuliaMath/HCubature.jl
// written in Julia by Steven G. Johnson served as a template.
// It comes with the following MIT "Expat" license:
//
// Copyright (c) 2017: Steven G. Johnson.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files
// (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

#ifndef STAN_MATH_PRIM_FUNCTOR_HCUBATURE_HPP
#define STAN_MATH_PRIM_FUNCTOR_HCUBATURE_HPP

#include <stan/math/prim/fun/choose.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/fmax.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <queue>
#include <tuple>
#include <iostream>

namespace stan {
namespace math {

namespace internal {

static constexpr std::array<double, 8> xd7{
    -9.9145537112081263920685469752598e-01,
    -9.4910791234275852452618968404809e-01,
    -8.6486442335976907278971278864098e-01,
    -7.415311855993944398638647732811e-01,
    -5.8608723546769113029414483825842e-01,
    -4.0584515137739716690660641207707e-01,
    -2.0778495500789846760068940377309e-01,
    0.0};

static constexpr std::array<double, 8> wd7{
    2.2935322010529224963732008059913e-02,
    6.3092092629978553290700663189093e-02,
    1.0479001032225018383987632254189e-01,
    1.4065325971552591874518959051021e-01,
    1.6900472663926790282658342659795e-01,
    1.9035057806478540991325640242055e-01,
    2.0443294007529889241416199923466e-01,
    2.0948214108472782801299917489173e-01};

static constexpr std::array<double, 4> gwd7{
    1.2948496616886969327061143267787e-01, 2.797053914892766679014677714229e-01,
    3.8183005050511894495036977548818e-01,
    4.1795918367346938775510204081658e-01};

/**
 * Get the [x]-th lexicographically ordered set of [p] elements in [dim]
 * output is in [c], and should be sizeof(int)*[p]
 * "Algorithm 515: Generation of a Vector from the Lexicographical Index";
 * Buckles, B. P., and Lybanon, M. ACM Transactions on Mathematical Software,
 * Vol. 3, No. 2, June 1977. User lucaroni from
 * https://stackoverflow.com/questions/561/how-to-use-combinations-of-sets-as-test-data#794
 *
 * @param[out] c output vector
 * @param dim dimension
 * @param p number of elements
 * @param x x-th lexicographically ordered set
 */
inline Eigen::VectorXi combination(Eigen::VectorXi& c, const int dim,
                                   const int p, const int x) {
  int r = 0;
  int k = 0;
  c[0] = 0;
  for (; k < x; k = k + r) {
    c[0]++;
    r = choose(dim - c[0], p - 1);
  }
  k = k - r;
  for (int i = 1; i < p - 1; i++) {
    c[i] = c[i - 1];
    for (; k < x; k = k + r) {
      c[i]++;
      r = choose(dim - c[i], p - (i + 1));
    }
    k = k - r;
  }
  if (p > 1) {
    c[p - 1] = c[p - 2] + x - k;
  } else {
    c[0] = x;
  }
  return std::move(c);
}

/**
 * Compute a vector [p] of all [dim]-component vectors
 * with [k] components equal to [lambda] and other components equal to zero.
 *
 * @param[in,out] p vector of vectors
 * @param k number of components equal to lambda
 * @param lambda scalar
 * @param dim dimension
 */
inline Eigen::MatrixXd combos(const int k, const double lambda, const int dim) {
  Eigen::VectorXi c(k);
  const auto choose_dimk = choose(dim, k);
  Eigen::MatrixXd p = Eigen::MatrixXd::Zero(dim, choose_dimk);
  for (size_t i = 0; i < choose_dimk; i++) {
    c = combination(c, dim, k, i + 1);
    for (size_t j = 0; j < k; j++) {
      p.coeffRef(c.coeff(j) - 1, i) = lambda;
    }
  }
  return p;
}

/**
 * Compute a vector [p] of all [dim]-component vectors
 * with [k] components equal to [±lambda] and other components equal to zero
 * (with all possible signs).
 *
 * @param[in,out] p vector of vectors
 * @param k number of components equal to lambda
 * @param lambda scalar
 * @param dim dimension
 */
inline Eigen::MatrixXd signcombos(const int k, const double lambda,
                                  const int dim) {
  Eigen::VectorXi c(k);
  const auto choose_dimk = choose(dim, k);
  Eigen::MatrixXd p = Eigen::MatrixXd::Zero(dim, choose_dimk * std::pow(2, k));
  int current_col = 0;
  for (int i = 1; i != choose_dimk + 1; i++) {
    c = combination(c, dim, k, i);
    std::vector<bool> index;
    for (int j = 0; j != std::pow(2, k); j++) {
      int prev_col = (j == 0) ? current_col : current_col - 1;
      p.col(current_col) = p.col(prev_col);

      if (index.size() == 0) {
        index.push_back(false);
        for (int h = 0; h != k; h++) {
          p.col(current_col)[c[h] - 1] = lambda;
        }
      } else {
        int first_zero = 0;
        while ((first_zero < index.size()) && index[first_zero]) {
          first_zero++;
        }
        if (first_zero == index.size()) {
          index.flip();
          for (int h = 0; h != index.size(); h++) {
            p.col(current_col)[c[h] - 1] *= -1;
          }
          index.push_back(true);
          p.col(current_col)[c[index.size() - 1] - 1] = -lambda;
        } else {
          for (int h = 0; h != first_zero + 1; h++) {
            if (index[h]) {
              index[h] = 0;
            } else {
              index[h] = 1;
            }
            p.col(current_col)[c[h] - 1] *= -1;
          }
        }
      }
      current_col += 1;
    }
  }
  return p;
}

/**
 * Compute the integral of the function to be integrated (integrand) from a to b
 * for one dimension.
 *
 * @tparam F type of the integrand
 * @tparam ParsPairT type of the pair of parameters for the integrand
 * @param integrand function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param pars_pair Pair of parameters for the integrand
 * @return numeric integral of the integrand and error
 */
template <typename F, typename ParsPairT>
std::pair<double, double> gauss_kronrod(const F& integrand, const double a,
                                        const double b,
                                        const ParsPairT& pars_pair) {
  Eigen::VectorXd c{{0.5 * (a + b)}};
  Eigen::VectorXd cp(1);
  Eigen::VectorXd cm(1);
  double delta = 0.5 * (b - a);
  double f0 = math::apply(
      [&integrand, &c](auto&&... args) { return integrand(c, args...); },
      pars_pair);

  double I = f0 * wd7[7];
  double Idash = f0 * gwd7[3];
  for (std::size_t i = 0; i != 7; i++) {
    double deltax = delta * xd7[i];
    cp[0] = c[0] + deltax;
    cm[0] = c[0] - deltax;
    double fx = math::apply(
        [&integrand, &cp](auto&&... args) { return integrand(cp, args...); },
        pars_pair);
    double temp = math::apply(
        [&integrand, &cm](auto&&... args) { return integrand(cm, args...); },
        pars_pair);
    fx += temp;
    I += fx * wd7[i];
    if (i % 2 == 1) {
      Idash += fx * gwd7[i / 2];
    }
  }
  double V = fabs(delta);
  I *= V;
  Idash *= V;
  return std::make_pair(I, fabs(I - Idash));
}

/**
 * Compute the points and weights corresponding to a [dim]-dimensional
 * Genz-Malik cubature rule
 *
 * @param[in,out] points points for the last 4 GenzMalik weights
 * @param[in,out] weights weights for the 5 terms in the GenzMalik rule
 * @param[in,out] weights_low_deg weights for the embedded lower-degree rule
 * @param dim dimension
 */
inline std::tuple<std::vector<Eigen::MatrixXd>, Eigen::VectorXd,
                  Eigen::VectorXd>
make_GenzMalik(const int dim) {
  std::vector<Eigen::MatrixXd> points(4);
  Eigen::VectorXd weights(5);
  Eigen::VectorXd weights_low_deg(4);
  double l4 = std::sqrt(9 * 1.0 / 10);
  double l2 = std::sqrt(9 * 1.0 / 70);
  double l3 = l4;
  double l5 = std::sqrt(9 * 1.0 / 19);

  double twopn = std::pow(2, dim);

  weights[0] = twopn * ((12824 - 9120 * dim + 400 * dim * dim) * 1.0 / 19683);
  weights[1] = twopn * (980.0 / 6561);
  weights[2] = twopn * ((1820 - 400 * dim) * 1.0 / 19683);
  weights[3] = twopn * (200.0 / 19683);
  weights[4] = 6859.0 / 19683;
  weights_low_deg[3] = twopn * (25.0 / 729);
  weights_low_deg[2] = twopn * ((265 - 100 * dim) * 1.0 / 1458);
  weights_low_deg[1] = twopn * (245.0 / 486);
  weights_low_deg[0] = twopn * ((729 - 950 * dim + 50 * dim * dim) * 1.0 / 729);
  points[0] = combos(1, l2, dim);
  points[1] = combos(1, l3, dim);
  points[2] = signcombos(2, l4, dim);
  points[3] = signcombos(dim, l5, dim);
  return std::make_tuple(std::move(points), std::move(weights),
                         std::move(weights_low_deg));
}

/**
 * Compute the integral of the function to be integrated (integrand) from a to b
 * for more than one dimensions.
 *
 * @tparam F type of the integrand
 * @tparam ParsTupleT type of the tuple of parameters for the integrand
 * @param[out] integrand function to be integrated
 * @param[in] points points for the last 4 GenzMalik weights
 * @param[in] weights weights for the 5 terms in the GenzMalik rule
 * @param[in] weights_low_deg weights for the embedded lower-degree rule
 * @param dim dimension of the multidimensional integral
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param pars_tuple Tuple of parameters for the integrand
 * @return numeric integral of the integrand, error, and suggested coordinate to
 * subdivide next
 */
template <typename F, typename ParsTupleT>
std::tuple<double, double, int> integrate_GenzMalik(
    const F& integrand, const std::vector<Eigen::MatrixXd>& points,
    const Eigen::VectorXd& weights, const Eigen::VectorXd& weights_low_deg,
    const int dim, const Eigen::VectorXd& a, const Eigen::VectorXd& b,
    const ParsTupleT& pars_tuple) {
  Eigen::VectorXd c = Eigen::VectorXd::Zero(dim);
  for (std::size_t i = 0; i != dim; i++) {
    if (a[i] == b[i]) {
      return std::make_tuple(0.0, 0.0, 0);
    }
    c[i] = (a[i] + b[i]) / 2;
  }
  Eigen::VectorXd deltac = ((b - a).array() / 2.0).matrix();
  double v = 1.0;
  for (std::size_t i = 0; i != dim; i++) {
    v *= deltac[i];
  }

  double f1 = math::apply(
      [&integrand, &c](auto&&... args) { return integrand(c, args...); },
      pars_tuple);
  double f2 = 0.0;
  double f3 = 0.0;
  double twelvef1 = 12 * f1;

  double maxdivdiff = 0.0;
  Eigen::VectorXd divdiff(dim);
  Eigen::VectorXd p2(dim);
  Eigen::VectorXd p3(dim);
  Eigen::VectorXd cc(dim);

  for (std::size_t i = 0; i != dim; i++) {
    for (std::size_t j = 0; j != dim; j++) {
      p2[j] = deltac[j] * points[0](j, i);
    }

    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p2[j];
    }
    double f2i = math::apply(
        [&integrand, &cc](auto&&... args) { return integrand(cc, args...); },
        pars_tuple);
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] - p2[j];
    }
    double temp = math::apply(
        [&integrand, &cc](auto&&... args) { return integrand(cc, args...); },
        pars_tuple);
    f2i += temp;

    for (std::size_t j = 0; j != dim; j++) {
      p3[j] = deltac[j] * points[1](j, i);
    }
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p3[j];
    }
    double f3i = math::apply(
        [&integrand, &cc](auto&&... args) { return integrand(cc, args...); },
        pars_tuple);
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] - p3[j];
    }
    temp = math::apply(
        [&integrand, &cc](auto&&... args) { return integrand(cc, args...); },
        pars_tuple);
    f3i += temp;
    f2 += f2i;
    f3 += f3i;
    divdiff[i] = fabs(f3i + twelvef1 - 7 * f2i);
  }
  Eigen::VectorXd p4(dim);
  double f4 = 0.0;
  for (std::size_t i = 0; i != points[2].cols(); i++) {
    for (std::size_t j = 0; j != dim; j++) {
      p4[j] = deltac[j] * points[2](j, i);
    }
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p4[j];
    }
    double temp = math::apply(
        [&integrand, &cc](auto&&... args) { return integrand(cc, args...); },
        pars_tuple);
    f4 += temp;
  }
  double f5 = 0.0;
  Eigen::VectorXd p5(dim);
  for (std::size_t i = 0; i != points[3].cols(); i++) {
    for (std::size_t j = 0; j != dim; j++) {
      p5[j] = deltac[j] * points[3](j, i);
    }

    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p5[j];
    }
    double temp = math::apply(
        [&integrand, &cc](auto&&... args) { return integrand(cc, args...); },
        pars_tuple);
    f5 += temp;
  }

  double I = v
             * (weights[0] * f1 + weights[1] * f2 + weights[2] * f3
                + weights[3] * f4 + weights[4] * f5);
  double Idash = v
                 * (weights_low_deg[0] * f1 + weights_low_deg[1] * f2
                    + weights_low_deg[2] * f3 + weights_low_deg[3] * f4);
  double E = fabs(I - Idash);

  int kdivide = 0;
  double deltaf = E / (std::pow(10, dim) * v);
  for (std::size_t i = 0; i != dim; i++) {
    double delta = divdiff[i] - maxdivdiff;
    if (delta > deltaf) {
      kdivide = i;
      maxdivdiff = divdiff[i];
    } else if ((fabs(delta) <= deltaf) && (deltac[i] > deltac[kdivide])) {
      kdivide = i;
    }
  }
  return std::make_tuple(I, E, kdivide);
}

/**
 * Compute the integral of the function to be integrated (integrand) from a to b
 * for more than one dimensions.
 *
 * @tparam Vec1 Type of vector 1
 * @tparam Vec2 Type of vector 2
 * @param a lower bounds of the integral
 * @param b upper bounds of the integral
 * @param I value of the integral
 * @param kdivide number of subdividing the integration volume
 */
struct Box {
  template <typename Vec1, typename Vec2>
  Box(Vec1&& a, Vec2&& b, double I, int kdivide)
      : a_(std::forward<Vec1>(a)),
        b_(std::forward<Vec2>(b)),
        I_(I),
        kdiv_(kdivide) {}
  Eigen::VectorXd a_;
  Eigen::VectorXd b_;
  double I_;
  int kdiv_;
};

}  // namespace internal

/**
 * Compute the dim-dimensional integral of the function \f$f\f$ from \f$a\f$ to
 \f$b\f$ within
 * specified relative and absolute tolerances or maximum number of evaluations.
 * \f$a\f$ and \f$b\f$ can be finite or infinite and should be given as vectors.
 *
 * The prototype for \f$f\f$ is:
 \verbatim
   template <typename T_x, typename T_p>
   stan::return_type_t<T_x, T_p> f(const T_x& x, const T_p& p) {
   using T_x_ref = stan::ref_type_t<T_x>;
   T_x_ref x_ref = x;
   stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
   my_params* pars = static_cast<my_params*>(p);
   type_1 var_1 = (pars->par_1);
   return ;
   }
 \endverbatim
 *
 * The parameters can be handed over to f as a struct:
 \verbatim
  struct my_params {
  type_1 par_1;
  type_2 par_2;
  };
 \endverbatim
 *
 * @tparam F Type of f
 * @tparam T_pars Type of paramete-struct
 *
 * @param integrand a functor with signature given above
 * @param pars parameters to be passed to f as a struct
 * @param dim dimension of the integral
 * @param a lower limit of integration as vector
 * @param b upper limit of integration as vector
 * @param max_eval maximal number of evaluations
 * @param reqAbsError absolute error
 * @param reqRelError relative error as vector
 * @param val correct value of integral
 *
 * @return The value of the dim-dimensional integral of \f$f\f$ from \f$a\f$ to
 \f$b\f$.
 * @throw std::domain_error no errors will be thrown.
 */
template <typename F, typename ParsTuple>
double hcubature(const F& integrand, const ParsTuple& pars, const int dim,
                 const Eigen::VectorXd& a, const Eigen::VectorXd& b,
                 const int max_eval, const double reqAbsError,
                 const double reqRelError) {
  const auto maxEval = max_eval <= 0 ? 1000000 : max_eval;
  double result;
  double err;
  int kdivide = 0;

  std::vector<Eigen::MatrixXd> p(4);
  Eigen::VectorXd w_five(5);
  Eigen::VectorXd wd_four(4);

  if (dim == 1) {
    std::tie(result, err)
        = internal::gauss_kronrod(integrand, a[0], b[0], pars);
  } else {
    std::tie(p, w_five, wd_four) = internal::make_GenzMalik(dim);
    std::tie(result, err, kdivide) = internal::integrate_GenzMalik(
        integrand, p, w_five, wd_four, dim, a, b, pars);
  }
  int numevals
      = (dim == 1) ? 15 : 1 + 4 * dim + 2 * dim * (dim - 1) + std::pow(2, dim);
  int evals_per_box = numevals;

  if ((err <= fmax(reqRelError * fabs(result), reqAbsError))
      || (numevals >= maxEval)) {
    return result;
  }
  numevals += 2 * evals_per_box;
  std::vector<internal::Box> ms;
  ms.reserve(numevals);
  ms.emplace_back(std::move(a), std::move(b), result, kdivide);
  auto get_largest_box_idx = [](auto&& box_vec) {
    auto max_it = std::max_element(box_vec.begin(), box_vec.end());
    return std::distance(box_vec.begin(), max_it);
  };
  std::vector<double> err_vec;
  err_vec.reserve(numevals);
  err_vec.push_back(err);
  while ((numevals < maxEval)
         && (err > max(reqRelError * fabs(result), reqAbsError))
         && std::isfinite(result)) {
    auto err_idx = get_largest_box_idx(err_vec);
    auto&& box = ms[err_idx];

    double w = (box.b_[box.kdiv_] - box.a_[box.kdiv_]) / 2;
    Eigen::VectorXd ma
        = Eigen::Map<const Eigen::VectorXd>(box.a_.data(), box.a_.size());

    ma[box.kdiv_] += w;
    Eigen::VectorXd mb
        = Eigen::Map<const Eigen::VectorXd>(box.b_.data(), box.b_.size());
    mb[box.kdiv_] -= w;

    double result_1;
    double result_2;
    double err_1;
    double err_2;
    double kdivide_1{0};
    double kdivide_2{0};

    if (dim == 1) {
      std::tie(result_1, err_1)
          = internal::gauss_kronrod(integrand, ma[0], box.b_[0], pars);
      std::tie(result_2, err_2)
          = internal::gauss_kronrod(integrand, box.a_[0], mb[0], pars);
    } else {
      std::tie(result_1, err_1, kdivide_1) = internal::integrate_GenzMalik(
          integrand, p, w_five, wd_four, dim, ma, box.b_, pars);
      std::tie(result_2, err_2, kdivide_2) = internal::integrate_GenzMalik(
          integrand, p, w_five, wd_four, dim, box.a_, mb, pars);
    }
    internal::Box box1(std::move(ma), std::move(box.b_), result_1, kdivide_1);
    internal::Box box2(std::move(box.a_), std::move(mb), result_2, kdivide_2);
    result += result_1 + result_2 - box.I_;
    err += err_1 + err_2 - err_vec[err_idx];
    ms[err_idx].I_ = 0;
    err_vec[err_idx] = 0;
    ms.push_back(std::move(box1));
    ms.push_back(std::move(box2));
    err_vec.push_back(err_1);
    err_vec.push_back(err_2);
    numevals += 2 * evals_per_box;
  }
  result = 0.0;
  err = 0.0;

  for (auto&& box : ms) {
    result += box.I_;
  }
  for (auto err_i : err_vec) {
    err += err_i;
  }
  return result;
}  // hcubature

}  // namespace math
}  // namespace stan
#endif
