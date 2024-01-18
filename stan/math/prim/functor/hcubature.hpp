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
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace stan {
namespace math {

namespace internal {

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
inline void combination(Eigen::Matrix<int, Eigen::Dynamic, 1>& c, const int dim,
                        const int p, const int x) {
  int r = 0;
  int k = 0;
  c[0] = 0;
  for (; k < x; r = choose(dim - c[0], p - 1), k = k + r) {
    c[0]++;
  }
  k = k - r;
  for (int i = 1; i < p - 1; i++) {
    c[i] = c[i - 1];
    for (; k < x; r = choose(dim - c[i], p - (i + 1)), k = k + r) {
      c[i]++;
    }
    k = k - r;
  }
  c[p - 1] = (p > 1) ? c[p - 2] + x - k : x;
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
template <typename Scalar>
inline Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> combos(
    const int k, const Scalar lambda, const int dim) {
  Eigen::Matrix<int, Eigen::Dynamic, 1> c(k);
  const auto choose_dimk = choose(dim, k);
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> p
      = Eigen::MatrixXd::Zero(dim, choose_dimk);
  for (size_t i = 0; i < choose_dimk; i++) {
    combination(c, dim, k, i + 1);
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
template <typename Scalar>
inline Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> signcombos(
    const int k, const Scalar lambda, const int dim) {
  Eigen::Matrix<int, Eigen::Dynamic, 1> c(k);
  const auto choose_dimk = choose(dim, k);
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> p
      = Eigen::MatrixXd::Zero(dim, choose_dimk * std::pow(2, k));
  int current_col = 0;
  const auto inner_iter_len = std::pow(2, k);
  std::vector<bool> index;
  index.reserve(inner_iter_len * (choose_dimk + 1));
  for (int i = 1; i != choose_dimk + 1; i++) {
    combination(c, dim, k, i);
    index.push_back(false);
    p(c.array() - 1.0, current_col).setConstant(lambda);
    current_col += 1;
    for (int j = 1; j != inner_iter_len; j++, current_col++) {
      p.col(current_col) = p.col(current_col - 1);
      int first_zero
          = std::distance(std::begin(index),
                          std::find(std::begin(index), std::end(index), false));
      const std::size_t index_size = index.size();
      if (first_zero == index_size) {
        index.flip();
        p(c.segment(0, index.size()).array() - 1, current_col).array() *= -1;
        index.push_back(true);
        p(c[index.size() - 1] - 1, current_col) = -lambda;
      } else {
        for (int h = 0; h != first_zero + 1; h++) {
          index[h].flip();
        }
        p(c.segment(0, first_zero + 1).array() - 1, current_col).array() *= -1;
      }
    }
    index.clear();
  }
  return p;
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

inline std::tuple<
    std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 4>,
    Eigen::Matrix<double, 5, 1>, Eigen::Matrix<double, 4, 1>>
make_GenzMalik(const int dim) {
  std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 4> points;
  Eigen::Matrix<double, 5, 1> weights;
  Eigen::Matrix<double, 4, 1> weights_low_deg;
  double twopn = std::pow(2, dim);
  weights[0]
      = twopn * ((12824.0 - 9120.0 * dim + 400.0 * dim * dim) * 1.0 / 19683.0);
  weights[1] = twopn * (980.0 / 6561.0);
  weights[2] = twopn * ((1820.0 - 400.0 * dim) * 1.0 / 19683.0);
  weights[3] = twopn * (200.0 / 19683.0);
  weights[4] = 6859.0 / 19683.0;
  weights_low_deg[0] = twopn * ((729 - 950 * dim + 50 * dim * dim) * 1.0 / 729);
  weights_low_deg[1] = twopn * (245.0 / 486);
  weights_low_deg[2] = twopn * ((265.0 - 100.0 * dim) * 1.0 / 1458.0);
  weights_low_deg[3] = twopn * (25.0 / 729.0);
  points[0] = combos(1, std::sqrt(9.0 * 1.0 / 70.0), dim);
  double l3 = std::sqrt(9.0 * 1.0 / 10.0);
  points[1] = combos(1, l3, dim);
  points[2] = signcombos(2, l3, dim);
  points[3] = signcombos(dim, std::sqrt(9.0 * 1.0 / 19.0), dim);
  return std::make_tuple(std::move(points), std::move(weights),
                         std::move(weights_low_deg));
}

/**
 * Compute the integral of the function to be integrated (integrand) from a to b
 * for more than one dimensions.
 *
 * @tparam F type of the integrand
 * @tparam GenzMalik type of -----
 * @tparam T_a type of lower limit of integration
 * @tparam T_b type of upper limit of integration
 * @tparam ParsTupleT type of the tuple of parameters for the integrand
 * @param[out] integrand function to be integrated
 * @param[in] genz_malik tuple of GenzMalik weights
 * @param dim dimension of the multidimensional integral
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param pars_tuple Tuple of parameters for the integrand
 * @return numeric integral of the integrand, error, and suggested coordinate to
 * subdivide next
 */
template <typename F, typename GenzMalik, typename T_a, typename T_b,
          typename ParsTupleT>
inline auto integrate_GenzMalik(const F& integrand, const GenzMalik& genz_malik,
                                const int dim,
                                const Eigen::Matrix<T_a, Eigen::Dynamic, 1>& a,
                                const Eigen::Matrix<T_b, Eigen::Dynamic, 1>& b,
                                const ParsTupleT& pars_tuple) {
  auto&& points = std::get<0>(genz_malik);
  auto&& weights = std::get<1>(genz_malik);
  auto&& weights_low_deg = std::get<2>(genz_malik);
  using delta_t = return_type_t<T_a, T_b>;
  Eigen::Matrix<delta_t, Eigen::Dynamic, 1> c(dim);
  Eigen::Matrix<delta_t, Eigen::Dynamic, 1> deltac(dim);
  using Scalar = return_type_t<ParsTupleT, delta_t>;
  for (auto i = 0; i != dim; i++) {
    if (a[i] == b[i]) {
      return std::make_tuple(Scalar(0.0), Scalar(0.0), 0);
    }
    c[i] = (a[i] + b[i]) / 2;
    deltac[i] = (b[i] - a[i]) / 2.0;
  }
  delta_t v = 1.0;
  for (std::size_t i = 0; i != dim; i++) {
    v *= deltac[i];
  }
  Eigen::Matrix<Scalar, 5, 1> f = Eigen::Matrix<Scalar, 5, 1>::Zero();
  f.coeffRef(0)
      = math::apply([](auto&& integrand, auto&& c,
                       auto&&... args) { return integrand(c, args...); },
                    pars_tuple, integrand, c);
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> divdiff(dim);
  for (auto i = 0; i != dim; i++) {
    auto p2 = deltac.cwiseProduct(points[0].col(i));
    auto f2i = math::apply(
        [](auto&& integrand, auto&& c, auto&& p2, auto&&... args) {
          return integrand(c + p2, args...) + integrand(c - p2, args...);
        },
        pars_tuple, integrand, c, p2);
    auto p3 = deltac.cwiseProduct(points[1].col(i));
    auto f3i = math::apply(
        [](auto&& integrand, auto&& c, auto&& p3, auto&&... args) {
          return integrand(c + p3, args...) + integrand(c - p3, args...);
        },
        pars_tuple, integrand, c, p3);
    f.coeffRef(1) += f2i;
    f.coeffRef(2) += f3i;
    divdiff[i] = fabs(f3i + 12.0 * f.coeff(0) - 7.0 * f2i);
  }
  for (auto i = 0; i != points[2].cols(); i++) {
    f.coeffRef(3) += math::apply(
        [](auto&& integrand, auto&& cc, auto&&... args) {
          return integrand(cc, args...);
        },
        pars_tuple, integrand, c + deltac.cwiseProduct(points[2].col(i)));
  }
  for (auto i = 0; i != points[3].cols(); i++) {
    f.coeffRef(4) += math::apply(
        [](auto&& integrand, auto&& cc, auto&&... args) {
          return integrand(cc, args...);
        },
        pars_tuple, integrand, c + deltac.cwiseProduct(points[3].col(i)));
  }

  Scalar I = v * weights.dot(f);
  Scalar Idash = v * weights_low_deg.dot(f.template head<4>());
  Scalar E = fabs(I - Idash);

  int kdivide = 0;
  Scalar deltaf = E / (std::pow(10, dim) * v);
  Scalar maxdivdiff = 0.0;
  for (auto i = 0; i != dim; i++) {
    Scalar delta = divdiff[i] - maxdivdiff;
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
 * @tparam T_a Type of return_type_t 1
 * @tparam T_b Type of return_type_t 2
 * @param a lower bounds of the integral
 * @param b upper bounds of the integral
 * @param I value of the integral
 * @param kdivide number of subdividing the integration volume
 */
template <typename T_a, typename T_b>
struct Box {
  template <typename Vec1, typename Vec2>
  Box(Vec1&& a, Vec2&& b, return_type_t<T_a, T_b> I, int kdivide)
      : a_(std::forward<Vec1>(a)),
        b_(std::forward<Vec2>(b)),
        I_(I),
        kdiv_(kdivide) {}
  Eigen::Matrix<T_a, Eigen::Dynamic, 1> a_;
  Eigen::Matrix<T_b, Eigen::Dynamic, 1> b_;
  return_type_t<T_a, T_b> I_;
  int kdiv_;
};

}  // namespace internal

/**
 * Compute the [dim]-dimensional integral of the function \f$f\f$ from \f$a\f$
 to \f$b\f$ within
 * specified relative and absolute tolerances or maximum number of evaluations.
 * \f$a\f$ and \f$b\f$ can be finite or infinite and should be given as vectors.
 *
 * @tparam F Type of f
 * @tparam T_a Type of lower limit of integration
 * @tparam T_b Type of upper limit of integration
 * @tparam ParsTuple Type of parameter-tuple
 * @tparam TAbsErr Type of absolute error
 * @tparam TRelErr Type of relative error
 *
 * @param integrand a functor with signature given above
 * @param pars parameters to be passed to f as a tuple
 * @param dim dimension of the integral
 * @param a lower limit of integration as vector
 * @param b upper limit of integration as vector
 * @param max_eval maximal number of evaluations
 * @param reqAbsError absolute error
 * @param reqRelError relative error as vector
 *
 * @return The value of the [dim]-dimensional integral of \f$f\f$ from \f$a\f$
 to \f$b\f$.
 * @throw std::domain_error no errors will be thrown.
 */
template <typename F, typename T_a, typename T_b, typename ParsTuple,
          typename TAbsErr, typename TRelErr>
inline auto hcubature(const F& integrand, const ParsTuple& pars, const int dim,
                      const Eigen::Matrix<T_a, Eigen::Dynamic, 1>& a,
                      const Eigen::Matrix<T_b, Eigen::Dynamic, 1>& b,
                      const int max_eval, const TAbsErr reqAbsError,
                      const TRelErr reqRelError) {
  using Scalar = return_type_t<ParsTuple, T_a, T_b, TAbsErr, TRelErr>;
  using eig_vec_a = Eigen::Matrix<T_a, Eigen::Dynamic, 1>;
  using eig_vec_b = Eigen::Matrix<T_b, Eigen::Dynamic, 1>;
  using namespace boost::math::quadrature;

  const int maxEval = max_eval <= 0 ? 1000000 : max_eval;
  Scalar result;
  Scalar err;
  auto kdivide = 0;
  std::tuple<
      std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 4>,
      Eigen::Matrix<double, 5, 1>, Eigen::Matrix<double, 4, 1>>
      genz_malik;

  auto gk_lambda = [&integrand, &pars](auto&& c) {
    return stan::math::apply(
        [](auto&& integrand, auto&& c, auto&&... args) {
          return integrand(c, args...);
        },
        pars, integrand, c);
  };

  if (dim == 1) {
    result = gauss_kronrod<Scalar, 15>::integrate(gk_lambda, a[0], b[0], 15, 0,
                                                  &err);

  } else {
    genz_malik = internal::make_GenzMalik(dim);
    std::tie(result, err, kdivide)
        = internal::integrate_GenzMalik(integrand, genz_malik, dim, a, b, pars);
  }
  auto numevals
      = (dim == 1) ? 15 : 1 + 4 * dim + 2 * dim * (dim - 1) + std::pow(2, dim);

  auto evals_per_box = numevals;
  Scalar error = err;
  Scalar val = result;

  if ((error <= fmax(reqRelError * fabs(val), reqAbsError))
      || (numevals >= maxEval)) {
    return result;
  }
  numevals += 2 * evals_per_box;
  using box_t = internal::Box<T_a, T_b>;
  std::vector<box_t> ms;
  ms.reserve(numevals);
  ms.emplace_back(a, b, result, kdivide);
  auto get_largest_box_idx = [](auto&& box_vec) {
    auto max_it = std::max_element(box_vec.begin(), box_vec.end());
    return std::distance(box_vec.begin(), max_it);
  };
  std::vector<Scalar> err_vec;
  err_vec.reserve(numevals);
  err_vec.push_back(err);
  while ((numevals < maxEval)
         && (error > fmax(reqRelError * fabs(val), reqAbsError))
         && fabs(val) < INFTY) {
    auto err_idx = get_largest_box_idx(err_vec);
    auto&& box = ms[err_idx];
    auto w = (box.b_[box.kdiv_] - box.a_[box.kdiv_]) / 2;
    eig_vec_a ma = Eigen::Map<const eig_vec_a>(box.a_.data(), box.a_.size());
    ma[box.kdiv_] += w;
    eig_vec_b mb = Eigen::Map<const eig_vec_b>(box.b_.data(), box.b_.size());
    mb[box.kdiv_] -= w;
    int kdivide_1{0};
    int kdivide_2{0};
    Scalar result_1;
    Scalar result_2;
    Scalar err_1;
    Scalar err_2;
    if (dim == 1) {
      result_1 = gauss_kronrod<Scalar, 15>::integrate(gk_lambda, ma[0],
                                                      box.b_[0], 15, 0, &err_1);
      result_2 = gauss_kronrod<Scalar, 15>::integrate(gk_lambda, box.a_[0],
                                                      mb[0], 15, 0, &err_2);
    } else {
      std::tie(result_1, err_1, kdivide_1) = internal::integrate_GenzMalik(
          integrand, genz_malik, dim, ma, box.b_, pars);
      std::tie(result_2, err_2, kdivide_2) = internal::integrate_GenzMalik(
          integrand, genz_malik, dim, box.a_, mb, pars);
    }
    box_t box1(std::move(ma), std::move(box.b_), result_1, kdivide_1);
    box_t box2(std::move(box.a_), std::move(mb), result_2, kdivide_2);
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
