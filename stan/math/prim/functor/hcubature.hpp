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
#include <queue>
#include <tuple>

namespace stan {
namespace math {

namespace internal {

static constexpr double xd7[8] = {-9.9145537112081263920685469752598e-01,
                                  -9.4910791234275852452618968404809e-01,
                                  -8.6486442335976907278971278864098e-01,
                                  -7.415311855993944398638647732811e-01,
                                  -5.8608723546769113029414483825842e-01,
                                  -4.0584515137739716690660641207707e-01,
                                  -2.0778495500789846760068940377309e-01,
                                  0.0};

static constexpr double wd7[8] = {2.2935322010529224963732008059913e-02,
                                  6.3092092629978553290700663189093e-02,
                                  1.0479001032225018383987632254189e-01,
                                  1.4065325971552591874518959051021e-01,
                                  1.6900472663926790282658342659795e-01,
                                  1.9035057806478540991325640242055e-01,
                                  2.0443294007529889241416199923466e-01,
                                  2.0948214108472782801299917489173e-01};

static constexpr double gwd7[4] = {1.2948496616886969327061143267787e-01,
                                   2.797053914892766679014677714229e-01,
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
 * @param c output vector
 * @param dim dimension
 * @param p number of elements
 * @param x x-th lexicographically ordered set
 */
inline void combination(std::vector<int>& c, const int& dim, const int& p,
                        const int& x) {
  size_t r, k = 0;
  for (std::size_t i = 0; i < p - 1; i++) {
    c[i] = (i != 0) ? c[i - 1] : 0;
    do {
      c[i]++;
      r = choose(dim - c[i], p - (i + 1));
      k = k + r;
    } while (k < x);
    k = k - r;
  }
  if (p > 1) {
    c[p - 1] = c[p - 2] + x - k;
  } else {
    c[0] = x;
  }
}

/**
 * Compute a vector [p] of all [dim]-component vectors
 * with [k] components equal to [lambda] and other components equal to zero.
 *
 * @param k number of components equal to lambda
 * @param lambda scalar
 * @param dim dimension
 * @param p vector of vectors
 */
inline void combos(const int& k, const double& lambda, const int& dim,
                   std::vector<std::vector<double>>& p) {
  std::vector<int> c(k);
  const auto choose_dimk = choose(dim, k);
  for (std::size_t i = 1; i != choose_dimk + 1; i++) {
    std::vector<double> temp(dim, 0.0);
    combination(c, dim, k, i);
    for (std::size_t j = 0; j != k; j++) {
      temp[c[j] - 1] = lambda;
    }
    p.push_back(temp);
  }
}

/**
 * Helper function for signcombos
 *
 * @param index helper vector
 * @param k number of components equal to lambda
 * @param lambda scalar
 * @param c ordered vector
 * @param temp helper vector
 */
inline void increment(std::vector<bool>& index, const int& k,
                      const double& lambda, const std::vector<int>& c,
                      std::vector<double>& temp) {
  if (index.size() == 0) {
    index.push_back(false);
    for (std::size_t j = 0; j != k; j++) {
      temp[c[j] - 1] = lambda;
    }
    return;
  }
  int first_zero = 0;
  while ((first_zero < index.size()) && index[first_zero]) {
    first_zero++;
  }
  if (first_zero == index.size()) {
    index.flip();
    for (std::size_t j = 0; j != index.size(); j++) {
      temp[c[j] - 1] *= -1;
    }
    index.push_back(true);
    temp[c[index.size() - 1] - 1] = -lambda;
  } else {
    for (std::size_t i = 0; i != first_zero + 1; i++) {
      if (index[i]) {
        index[i] = 0;
      } else {
        index[i] = 1;
      }
      temp[c[i] - 1] *= -1;
    }
  }
}

/**
 * Compute a vector [p] of all [dim]-component vectors
 * with [k] components equal to [±lambda] and other components equal to zero
 * (with all possible signs).
 *
 * @param k number of components equal to lambda
 * @param lambda scalar
 * @param dim dimension
 * @param p vector of vectors
 */
inline void signcombos(const int& k, const double& lambda, const int& dim,
                       std::vector<std::vector<double>>& p) {
  std::vector<int> c(k);
  const auto choose_dimk = choose(dim, k);
  for (std::size_t i = 1; i != choose_dimk + 1; i++) {
    std::vector<double> temp(dim, 0.0);
    combination(c, dim, k, i);
    std::vector<bool> index;
    index.clear();
    for (std::size_t j = 0; j != std::pow(2, k); j++) {
      increment(index, k, lambda, c, temp);
      p.push_back(temp);
    }
  }
}

/**
 * Compute the integral of the function to be integrated (integrand) from a to b
 * for one dimension.
 *
 * @tparam F type of the integrand
 * @tparam T_pars type of the parameters for the integrand
 * @param integrand function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param pars parameters for the integrand
 * @return numeric integral of the integrand and error
 */
template <typename F, typename T_pars>
std::tuple<double, double> gauss_kronrod(const F& integrand, const double& a,
                                         const double& b, T_pars& pars) {
  std::vector<double> c(1, 0);
  std::vector<double> cp(1, 0);
  std::vector<double> cm(1, 0);
  c[0] = 0.5 * (a + b);
  double delta = 0.5 * (b - a);
  double f0 = integrand(c, pars);
  double I = f0 * wd7[7];
  double Idash = f0 * gwd7[3];
  for (std::size_t i = 0; i != 7; i++) {
    double deltax = delta * xd7[i];
    cp[0] = c[0] + deltax;
    cm[0] = c[0] - deltax;
    double fx = integrand(cp, pars);
    double temp = integrand(cm, pars);
    fx += temp;
    I += fx * wd7[i];
    if (i % 2 == 1) {
      Idash += fx * gwd7[i / 2];
    }
  }
  double V = fabs(delta);
  I *= V;
  Idash *= V;
  return std::make_tuple(I, fabs(I - Idash));
}

/**
 * Compute the points and weights corresponding to a [dim]-dimensional
 * Genz-Malik cubature rule
 *
 * @param dim dimension
 * @param p points for the last 4 GenzMalik weights
 * @param w weights for the 5 terms in the GenzMalik rule
 * @param wd weights for the embedded lower-degree rule
 */
inline void make_GenzMalik(const int& dim,
                           std::vector<std::vector<std::vector<double>>>& p,
                           std::vector<double>& w, std::vector<double>& wd) {
  double l4 = std::sqrt(9 * 1.0 / 10);
  double l2 = std::sqrt(9 * 1.0 / 70);
  double l3 = l4;
  double l5 = std::sqrt(9 * 1.0 / 19);

  double twopn = std::pow(2, dim);

  w[0] = twopn * ((12824 - 9120 * dim + 400 * dim * dim) * 1.0 / 19683);
  w[1] = twopn * (980.0 / 6561);
  w[2] = twopn * ((1820 - 400 * dim) * 1.0 / 19683);
  w[3] = twopn * (200.0 / 19683);
  w[4] = 6859.0 / 19683;
  wd[3] = twopn * (25.0 / 729);
  wd[2] = twopn * ((265 - 100 * dim) * 1.0 / 1458);
  wd[1] = twopn * (245.0 / 486);
  wd[0] = twopn * ((729 - 950 * dim + 50 * dim * dim) * 1.0 / 729);

  combos(1, l2, dim, p[0]);
  combos(1, l3, dim, p[1]);
  signcombos(2, l4, dim, p[2]);
  signcombos(dim, l5, dim, p[3]);
}

/**
 * Compute the integral of the function to be integrated (integrand) from a to b
 * for more than one dimensions.
 *
 * @tparam F type of the integrand
 * @tparam T_pars type of the parameters for the integrand
 * @param integrand function to be integrated
 * @param p
 * @param w
 * @param wd
 * @param dim dimension of the multidimensional integral
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param pars parameters for the integrand
 * @return numeric integral of the integrand, error, and suggested coordinate to
 * subdivide next
 */
template <typename F, typename T_pars>
std::tuple<double, double, int> integrate_GenzMalik(
    const F& integrand, std::vector<std::vector<std::vector<double>>>& p,
    std::vector<double>& w, std::vector<double>& wd, const int& dim,
    const std::vector<double>& a, const std::vector<double>& b, T_pars& pars) {
  std::vector<double> c(dim, 0);
  std::vector<double> deltac(dim);

  for (std::size_t i = 0; i != dim; i++) {
    c[i] = (a[i] + b[i]) / 2;
  }

  for (std::size_t i = 0; i != dim; i++) {
    deltac[i] = fabs(b[i] - a[i]) / 2;
  }
  double v = 1.0;
  for (std::size_t i = 0; i != dim; i++) {
    v *= deltac[i];
  }

  if (v == 0.0) {
    return std::make_tuple(0.0, 0.0, 0);
  }

  double f1 = integrand(c, pars);
  double f2 = 0.0;
  double f3 = 0.0;
  double twelvef1 = 12 * f1;

  double maxdivdiff = 0.0;
  std::vector<double> divdiff(dim);
  std::vector<double> p2(dim);
  std::vector<double> p3(dim);
  std::vector<double> cc(dim, 0);

  for (std::size_t i = 0; i != dim; i++) {
    for (std::size_t j = 0; j != dim; j++) {
      p2[j] = deltac[j] * p[0][i][j];
    }

    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p2[j];
    }
    double f2i = integrand(cc, pars);
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] - p2[j];
    }
    double temp = integrand(cc, pars);
    f2i += temp;

    for (std::size_t j = 0; j != dim; j++) {
      p3[j] = deltac[j] * p[1][i][j];
    }
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p3[j];
    }
    double f3i = integrand(cc, pars);
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] - p3[j];
    }
    temp = integrand(cc, pars);
    f3i += temp;
    f2 += f2i;
    f3 += f3i;
    divdiff[i] = fabs(f3i + twelvef1 - 7 * f2i);
  }
  std::vector<double> p4(dim);
  double f4 = 0.0;
  for (std::size_t i = 0; i != p[2].size(); i++) {
    for (std::size_t j = 0; j != dim; j++) {
      p4[j] = deltac[j] * p[2][i][j];
    }
    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p4[j];
    }
    double temp = integrand(cc, pars);
    f4 += temp;
  }
  double f5 = 0.0;
  std::vector<double> p5(dim);
  for (std::size_t i = 0; i != p[3].size(); i++) {
    for (std::size_t j = 0; j != dim; j++) {
      p5[j] = deltac[j] * p[3][i][j];
    }

    for (std::size_t j = 0; j != dim; j++) {
      cc[j] = c[j] + p5[j];
    }
    double temp = integrand(cc, pars);
    f5 += temp;
  }

  double I = v * (w[0] * f1 + w[1] * f2 + w[2] * f3 + w[3] * f4 + w[4] * f5);
  double Idash = v * (wd[0] * f1 + wd[1] * f2 + wd[2] * f3 + wd[3] * f4);
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

struct Box {
  Box(const std::vector<double>& a, const std::vector<double>& b, double I,
      double err, int kdivide)
      : a(a), b(b), I(I), E(err), kdiv(kdivide) {}
  bool operator<(const Box& box) const { return E < box.E; }
  std::vector<double> a;
  std::vector<double> b;
  double I;
  double E;
  int kdiv;
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
 * @param f a functor with signature given above
 * @param dim dimension of the integral
 * @param a lower limit of integration as vector
 * @param b upper limit of integration as vector
 * @param maxEval maximal number of evaluations
 * @param reqAbsError absolute error
 * @param reqRelError relative error as vector
 * @param pars parameters to be passed to f as a struct
 * @param val correct value of integral
 *
 * @return The value of the dim-dimensional integral of \f$f\f$ from \f$a\f$ to
 \f$b\f$.
 * @throw std::domain_error no errors will be thrown.
 */
template <typename F, typename T_pars>
double hcubature(const F& integrand, const T_pars& pars, const int& dim,
                 const std::vector<double>& a, const std::vector<double>& b,
                 int& maxEval, const double& reqAbsError,
                 const double& reqRelError) {
  if (maxEval <= 0) {
    maxEval = 1000000;
  }

  double result, err;
  int kdivide = 0;

  std::vector<std::vector<std::vector<double>>> p(4);
  std::vector<double> w_five(5);
  std::vector<double> wd_four(4);

  if (dim == 1) {
    std::tie(result, err)
        = internal::gauss_kronrod(integrand, a[0], b[0], pars);
  } else {
    internal::make_GenzMalik(dim, p, w_five, wd_four);
    std::tie(result, err, kdivide) = internal::integrate_GenzMalik(
        integrand, p, w_five, wd_four, dim, a, b, pars);
  }
  int numevals
      = (dim == 1) ? 15 : 1 + 4 * dim + 2 * dim * (dim - 1) + std::pow(2, dim);
  int evals_per_box = numevals;
  int kdiv = kdivide;
  double error = err;
  double val = result;

  if ((error <= fmax(reqRelError * fabs(val), reqAbsError))
      || (numevals >= maxEval)) {
    return val;
  }
  std::priority_queue<internal::Box> ms;
  internal::Box box(a, b, result, err, kdivide);
  ms.push(box);

  numevals += 2 * evals_per_box;
  while ((numevals < maxEval)
         && (error > max(reqRelError * fabs(val), reqAbsError))
         && std::isfinite(val)) {
    internal::Box box = ms.top();
    ms.pop();

    double w = (box.b[box.kdiv] - box.a[box.kdiv]) / 2;
    std::vector<double> ma(box.a);

    ma[box.kdiv] += w;
    std::vector<double> mb(box.b);
    mb[box.kdiv] -= w;

    double result_1, result_2, err_1, err_2, kdivide_1, kdivide_2;

    if (dim == 1) {
      std::tie(result_1, err_1)
          = internal::gauss_kronrod(integrand, ma[0], box.b[0], pars);
      std::tie(result_2, err_2)
          = internal::gauss_kronrod(integrand, box.a[0], mb[0], pars);
    } else {
      std::tie(result_1, err_1, kdivide_1) = internal::integrate_GenzMalik(
          integrand, p, w_five, wd_four, dim, ma, box.b, pars);
      std::tie(result_2, err_2, kdivide_2) = internal::integrate_GenzMalik(
          integrand, p, w_five, wd_four, dim, box.a, mb, pars);
    }
    internal::Box box1(ma, box.b, result_1, err_1, kdivide_1);
    ms.push(box1);
    internal::Box box2(box.a, mb, result_2, err_2, kdivide_2);
    ms.push(box2);
    val += box1.I + box2.I - box.I;
    error += box1.E + box2.E - box.E;
    numevals += 2 * evals_per_box;
  }
  val = 0.0;
  error = 0.0;

  for (; !ms.empty(); ms.pop()) {
    internal::Box box = ms.top();
    val += box.I;
    error += box.E;
  }
  return val;
}  // hcubature

}  // namespace math
}  // namespace stan
#endif
