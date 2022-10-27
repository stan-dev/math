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

#include <set>
#include <stan/math/prim/fun/choose.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/fmax.hpp>
#include <stan/math/prim/fun/max.hpp>

namespace stan {
namespace math {

namespace internal {

// tools
static const double xd7[8] = {-9.9145537112081263920685469752598e-01,
                              -9.4910791234275852452618968404809e-01,
                              -8.6486442335976907278971278864098e-01,
                              -7.415311855993944398638647732811e-01,
                              -5.8608723546769113029414483825842e-01,
                              -4.0584515137739716690660641207707e-01,
                              -2.0778495500789846760068940377309e-01,
                              0.0};

static const double wd7[8] = {2.2935322010529224963732008059913e-02,
                              6.3092092629978553290700663189093e-02,
                              1.0479001032225018383987632254189e-01,
                              1.4065325971552591874518959051021e-01,
                              1.6900472663926790282658342659795e-01,
                              1.9035057806478540991325640242055e-01,
                              2.0443294007529889241416199923466e-01,
                              2.0948214108472782801299917489173e-01};

static const double gwd7[4] = {1.2948496616886969327061143267787e-01,
                               2.797053914892766679014677714229e-01,
                               3.8183005050511894495036977548818e-01,
                               4.1795918367346938775510204081658e-01};

struct one_d {
  double result;
  double err;
  int kdivide = 0;
};

struct GenzMalik {
  std::vector<double*> p[4];
  double w[5];
  double wd[4];
};

template <typename T_c, typename T_n, typename T_p, typename T_x>
void combination(T_c& c, const T_n& n, const T_p& p, const T_x& x) {
  size_t i, r, k = 0;
  for (i = 0; i < p - 1; i++) {
    c[i] = (i != 0) ? c[i - 1] : 0;
    do {
      c[i]++;
      r = choose(n - c[i], p - (i + 1));
      k = k + r;
    } while (k < x);
    k = k - r;
  }
  if (p > 1)
    c[p - 1] = c[p - 2] + x - k;
  else
    c[0] = x;
}

template <typename T_k, typename T_lambda, typename T_n, typename T_p>
void combos(const T_k& k, const T_lambda& lambda, const T_n& n, T_p& p) {
  T_k* c = reinterpret_cast<T_k*>(malloc(k * sizeof(T_k)));
  for (size_t i = 1; i != choose(n, k) + 1; i++) {
    T_lambda* temp = reinterpret_cast<T_lambda*>(calloc(n, sizeof(T_lambda)));
    combination(c, n, k, i);
    for (size_t j = 0; j != k; j++)
      temp[c[j] - 1] = lambda;
    p.push_back(temp);
  }
  free(c);
}

template <typename T_index, typename T_k, typename T_lambda, typename T_n,
          typename T_c, typename T_temp>
void increment(T_index& index, const T_k& k, const T_lambda& lambda,
               const T_n& n, const T_c& c, T_temp& temp) {
  // temp size n, all elements initially zero
  if (index.size() == 0) {
    index.push_back(false);
    for (size_t j = 0; j != k; j++)
      temp[c[j] - 1] = lambda;
    return;
  }
  T_k first_zero = 0;
  while ((first_zero < index.size()) && index[first_zero])
    first_zero++;
  if (first_zero == index.size()) {
    index.flip();
    for (size_t j = 0; j != index.size(); j++)
      temp[c[j] - 1] *= -1;
    index.push_back(true);
    temp[c[index.size() - 1] - 1] = -lambda;
  } else {
    for (size_t i = 0; i != first_zero + 1; i++) {
      //           index[i] = !index[i];
      if (index[i])
        index[i] = 0;
      else
        index[i] = 1;
      temp[c[i] - 1] *= -1;
    }
  }
}

template <typename T_k, typename T_lambda, typename T_n, typename T_p>
void signcombos(const T_k& k, const T_lambda& lambda, const T_n& n, T_p& p) {
  T_k* c = reinterpret_cast<T_k*>(malloc(k * sizeof(T_k)));
  for (size_t i = 1; i != choose(n, k) + 1; i++) {
    T_lambda* temp = reinterpret_cast<T_lambda*>(calloc(n, sizeof(T_lambda)));
    combination(c, n, k, i);
    std::vector<bool> index;
    index.clear();
    for (size_t j = 0; j != std::pow(2, k); j++) {
      increment(index, k, lambda, n, c, temp);
      T_lambda* next
          = reinterpret_cast<T_lambda*>(malloc(n * sizeof(T_lambda)));
      memcpy(next, temp, n * sizeof(T_lambda));
      p.push_back(next);
    }
    free(temp);
  }
  free(c);
}

template <typename F, typename T_a, typename T_b, typename T_oned,
          typename T_pars>
void gauss_kronrod(const F& integrand, const T_a& a, const T_b& b, T_oned& out,
                   T_pars& pars) {
  using T_return_type = return_type_t<T_a, T_b, T_oned, T_pars>;
  std::vector<T_return_type> c(1, 0);
  std::vector<T_return_type> cp(1, 0);
  std::vector<T_return_type> cm(1, 0);
  c[0] = 0.5 * (a + b);
  T_return_type delta = 0.5 * (b - a);
  T_return_type f0 = integrand(c, pars);
  T_return_type I = f0 * wd7[7];
  T_return_type Idash = f0 * gwd7[3];
  for (size_t i = 0; i != 7; i++) {
    T_return_type deltax = delta * xd7[i];
    cp[0] = c[0] + deltax;
    cm[0] = c[0] - deltax;
    T_return_type fx = integrand(cp, pars);
    T_return_type temp = integrand(cm, pars);
    fx += temp;
    I += fx * wd7[i];
    if (i % 2 == 1)
      Idash += fx * gwd7[i / 2];
  }
  T_return_type V = fabs(delta);
  I *= V;
  Idash *= V;
  out.result = I;
  out.err = fabs(I - Idash);
}

template <typename T_n, typename T_GenzMalik>
void make_GenzMalik(const T_n& n, T_GenzMalik& g) {
  using T_return_type = return_type_t<T_n, T_GenzMalik>;
  T_return_type l4 = std::sqrt(9 * 1.0 / 10);
  T_return_type l2 = std::sqrt(9 * 1.0 / 70);
  T_return_type l3 = l4;
  T_return_type l5 = std::sqrt(9 * 1.0 / 19);

  T_return_type twopn = std::pow(2, n);

  g.w[0] = twopn * ((12824 - 9120 * n + 400 * n * n) * 1.0 / 19683);
  g.w[1] = twopn * (980.0 / 6561);
  g.w[2] = twopn * ((1820 - 400 * n) * 1.0 / 19683);
  g.w[3] = twopn * (200.0 / 19683);
  g.w[4] = 6859.0 / 19683;
  g.wd[3] = twopn * (25.0 / 729);
  g.wd[2] = twopn * ((265 - 100 * n) * 1.0 / 1458);
  g.wd[1] = twopn * (245.0 / 486);
  g.wd[0] = twopn * ((729 - 950 * n + 50 * n * n) * 1.0 / 729);

  combos(1, l2, n, g.p[0]);
  combos(1, l3, n, g.p[1]);
  signcombos(2, l4, n, g.p[2]);
  signcombos(n, l5, n, g.p[3]);
}

template <typename T_GenzMalik>
void clean_GenzMalik(T_GenzMalik& g) {
  for (size_t j = 0; j != 4; j++)
    for (size_t i = 0; i != g.p[j].size(); i++)
      if (g.p[j][i])
        free(g.p[j][i]);
}

template <typename F, typename T_GenzMalik, typename T_n, typename T_a,
          typename T_b, typename T_oned, typename T_pars>
void integrate_GenzMalik(const F& integrand, T_GenzMalik& g, const T_n& n,
                         const T_a& a, const T_b& b, T_oned& out,
                         T_pars& pars) {
  using T_return_type = return_type_t<T_n, T_a, T_b>;
  std::vector<T_return_type> c(n, 0);
  double* deltac = reinterpret_cast<double*>(malloc(n * sizeof(double)));

  for (size_t i = 0; i != n; i++)
    c[i] = (a[i] + b[i]) / 2;

  for (size_t i = 0; i != n; i++)
    deltac[i] = fabs(b[i] - a[i]) / 2;
  T_return_type v = 1.0;
  for (size_t i = 0; i != n; i++)
    v *= deltac[i];

  if (v == 0.0) {
    out.err = 0.0;
    out.result = 0.0;
    out.kdivide = 0;
    free(deltac);
    return;
  }

  T_return_type f1 = integrand(c, pars);
  T_return_type f2 = 0.0;
  T_return_type f3 = 0.0;
  T_return_type twelvef1 = 12 * f1;

  T_return_type maxdivdiff = 0.0;
  T_return_type* divdiff
      = reinterpret_cast<T_return_type*>(malloc(n * sizeof(T_return_type)));
  T_return_type* p2
      = reinterpret_cast<T_return_type*>(malloc(n * sizeof(T_return_type)));
  T_return_type* p3
      = reinterpret_cast<T_return_type*>(malloc(n * sizeof(T_return_type)));
  std::vector<T_return_type> cc(n, 0);

  for (size_t i = 0; i != n; i++) {
    for (size_t j = 0; j != n; j++)
      p2[j] = deltac[j] * g.p[0][i][j];

    for (size_t j = 0; j != n; j++)
      cc[j] = c[j] + p2[j];
    T_return_type f2i = integrand(cc, pars);
    for (size_t j = 0; j != n; j++)
      cc[j] = c[j] - p2[j];
    T_return_type temp = integrand(cc, pars);
    f2i += temp;

    for (size_t j = 0; j != n; j++)
      p3[j] = deltac[j] * g.p[1][i][j];
    for (size_t j = 0; j != n; j++)
      cc[j] = c[j] + p3[j];
    T_return_type f3i = integrand(cc, pars);
    for (size_t j = 0; j != n; j++)
      cc[j] = c[j] - p3[j];
    temp = integrand(cc, pars);
    f3i += temp;
    f2 += f2i;
    f3 += f3i;
    divdiff[i] = fabs(f3i + twelvef1 - 7 * f2i);
  }
  free(p2);
  free(p3);
  T_return_type* p4
      = reinterpret_cast<T_return_type*>(malloc(n * sizeof(T_return_type)));
  T_return_type f4 = 0.0;
  for (size_t i = 0; i != g.p[2].size(); i++) {
    for (size_t j = 0; j != n; j++)
      p4[j] = deltac[j] * g.p[2][i][j];
    for (size_t j = 0; j != n; j++)
      cc[j] = c[j] + p4[j];
    T_return_type temp = integrand(cc, pars);
    f4 += temp;
  }
  free(p4);
  T_return_type f5 = 0.0;
  T_return_type* p5
      = reinterpret_cast<T_return_type*>(malloc(n * sizeof(T_return_type)));
  for (size_t i = 0; i != g.p[3].size(); i++) {
    for (size_t j = 0; j != n; j++)
      p5[j] = deltac[j] * g.p[3][i][j];

    for (size_t j = 0; j != n; j++)
      cc[j] = c[j] + p5[j];
    T_return_type temp = integrand(cc, pars);
    f5 += temp;
  }
  free(p5);
  T_return_type I
      = v
        * (g.w[0] * f1 + g.w[1] * f2 + g.w[2] * f3 + g.w[3] * f4 + g.w[4] * f5);
  T_return_type Idash
      = v * (g.wd[0] * f1 + g.wd[1] * f2 + g.wd[2] * f3 + g.wd[3] * f4);
  T_return_type E = fabs(I - Idash);

  int kdivide = 0;
  T_return_type deltaf = E / (std::pow(10, n) * v);
  for (size_t i = 0; i != n; i++) {
    T_return_type delta = divdiff[i] - maxdivdiff;
    if (delta > deltaf) {
      kdivide = i;
      maxdivdiff = divdiff[i];
    } else if ((fabs(delta) <= deltaf) && (deltac[i] > deltac[kdivide])) {
      kdivide = i;
    }
  }
  out.result = I;
  out.err = E;
  out.kdivide = kdivide;
  free(deltac);
  free(divdiff);
}

class Box {
 public:
  Box(std::vector<double> a, std::vector<double> b, double& I, double& err,
      int& kdivide)
      : a(a), b(b), I(I), E(err), kdiv(kdivide) {}
  bool operator<(const Box& box) const { return E > box.E; }
  std::vector<double> a;
  std::vector<double> b;
  double I;
  double E;
  int kdiv;
};

inline class Box make_box(int n, std::vector<double> a, std::vector<double> b,
                          one_d out) {
  std::vector<double> ac(a);
  std::vector<double> bc(b);
  Box box(ac, bc, out.result, out.err, out.kdivide);
  return box;
}

}  // namespace internal

/**
 * Compute the n-dimensional integral of the function \f$f\f$ from \f$a\f$ to
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
 * @tparam T_dim Type of dimension
 * @tparam T_a Type of a
 * @tparam T_b Type of b
 * @tparam T_maxEval Type of maximum number of evaluations
 * @tparam T_reqAbsError Type of absolute error
 * @tparam T_reqRelError Type of relative error
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
 * @return The value of the n-dimensional integral of \f$f\f$ from \f$a\f$ to
 \f$b\f$.
 * @throw std::domain_error no errors will be thrown.
 */

// hcubature
template <typename F, typename T_pars, typename T_dim, typename T_a, typename T_b,
          typename T_maxEval, typename T_reqAbsError, typename T_reqRelError>
return_type_t<T_dim, T_a, T_b, T_maxEval, T_reqAbsError, T_reqRelError> hcubature(
    const F& integrand, const T_pars& pars, const T_dim& dim, const T_a& a,
    const T_b& b, const T_maxEval& maxEval, const T_reqAbsError& reqAbsError,
    const T_reqRelError& reqRelError) {
  using T_return_type
      = return_type_t<T_dim, T_a, T_b, T_maxEval, T_reqAbsError, T_reqRelError>;

  using T_a_ref = ref_type_t<T_a>;
  using T_b_ref = ref_type_t<T_b>;

  T_a_ref a_ref = a;
  T_b_ref b_ref = b;

  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_b_ref> b_vec(b_ref);

  internal::one_d out;
  internal::GenzMalik g;

  if (dim == 1) {
    internal::gauss_kronrod(integrand, a_vec.val(0), b_vec.val(0), out, pars);
  } else {
    internal::make_GenzMalik(dim, g);
    internal::integrate_GenzMalik(integrand, g, dim, a, b, out, pars);
  }
  T_return_type numevals
      = (dim == 1) ? 15 : 1 + 4 * dim + 2 * dim * (dim - 1) + std::pow(2, dim);
  T_return_type evals_per_box = numevals;
  T_return_type kdiv = out.kdivide;
  T_return_type err = out.err;
  T_return_type val = out.result;
  // convergence test
  if ((err <= fmax(reqRelError * fabs(val), reqAbsError))
      || ((maxEval != 0) && (numevals >= maxEval))) {
    internal::clean_GenzMalik(g);
    return val;
  }

  std::multiset<internal::Box> ms;
  ms.insert(internal::make_box(dim, a, b, out));

  while (true) {
    std::multiset<internal::Box>::iterator it;
    it = ms.begin();
    internal::Box box = *it;
    ms.erase(it);
    // split along dimension kdiv
    T_return_type w = (box.b[box.kdiv] - box.a[box.kdiv]) / 2;
    std::vector<double> ma(box.a);

    ma[box.kdiv] += w;
    std::vector<double> mb(box.b);
    mb[box.kdiv] -= w;

    if (dim == 1) {
      internal::gauss_kronrod(integrand, ma[0], box.b[0], out, pars);
    } else {
      internal::integrate_GenzMalik(integrand, g, dim, ma, box.b, out, pars);
    }
    internal::Box box1 = make_box(dim, ma, box.b, out);
    ms.insert(box1);

    if (dim == 1) {
      internal::gauss_kronrod(integrand, box.a[0], mb[0], out, pars);
    } else {
      internal::integrate_GenzMalik(integrand, g, dim, box.a, mb, out, pars);
    }
    internal::Box box2 = make_box(dim, box.a, mb, out);
    ms.insert(box2);
    val += box1.I + box2.I - box.I;
    err += box1.E + box2.E - box.E;
    numevals += 2 * evals_per_box;

    if (((err <= max(reqRelError * fabs(val), reqAbsError))
         || ((maxEval != 0) && (numevals >= maxEval)))
        || !(std::isfinite(val))) {
      break;
    }
  }
  val = 0.0;
  err = 0.0;

  for (std::multiset<internal::Box>::reverse_iterator rit = ms.rbegin();
       rit != ms.rend(); rit++) {
    val += (*rit).I;
    err += (*rit).E;
  }

  internal::clean_GenzMalik(g);
  return val;
}  // hcubature

}  // namespace math
}  // namespace stan
#endif
