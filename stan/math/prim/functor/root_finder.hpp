#ifndef STAN_MATH_PRIM_FUNCTOR_ROOT_FINDER_HPP
#define STAN_MATH_PRIM_FUNCTOR_ROOT_FINDER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_bounded.hpp>
#include <boost/math/tools/roots.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {


// TODO: Can just have one signature with another function for constructing tuple

// When caller provides functions for value and derivative and second derivative
template <typename F, typename FDiv, typename FDivDiv, typename GuessScalar, typename MinScalar, typename MaxScalar, typename... Types,
  require_all_st_arithmetic<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr>
auto root_finder_tol(F&& f, FDiv&& f_div, FDivDiv&& f_div_div, GuessScalar guess, MinScalar min, MaxScalar max, int digits,
                 std::uintmax_t& max_iter, Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  return_type_t<GuessScalar> ret = 0;
  auto f_plus_div = [&f, &f_div, &f_div_div, &args...](auto&& g) {
    return std::make_tuple(f(g, args...), f_div(g, args...), f_div_div(g, args...));
  };
  try {
    ret = boost::math::tools::halley_iterate(f_plus_div,
       guess, min, max, digits, max_iter);
  } catch (const std::exception& e) {
    throw e;
  }
  return ret;
}
/*
// When caller provides functions for value and derivative
template <typename F, typename FDiv, typename GuessScalar, typename MinScalar, typename MaxScalar, typename... Types,
  require_all_st_arithmetic<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr>
double root_finder_tol(F&& f, FDiv&& f_div, GuessScalar guess, MinScalar min, MaxScalar max, int digits,
                 std::uintmax_t& max_iter, Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  double ret = 0;
  auto f_plus_div = [&f, &f_div, &args...](auto&& g) {
    return std::make_tuple(f(g, args...), f_div(g, args...));
  };
  try {
    ret = boost::math::tools::halley_iterate(f_plus_div,
       guess, min, max, digits, max_iter);
  } catch (const std::exception& e) {
    throw e;
  }
  return ret;
}

// When supplied function returns pair
template <typename F, typename GuessScalar, typename MinScalar, typename MaxScalar, typename... Types,
  require_all_st_arithmetic<GuessScalar, MinScalar, MaxScalar, Types...>* = nullptr>
double root_finder_tol(F&& f, GuessScalar guess, MinScalar min, MaxScalar max, int digits,
                 std::uintmax_t& max_iter, Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  double ret = 0;
  auto f_plus_args = [&f, &args...](auto&& g) { return f(g, args...);};
  try {
    ret = boost::math::tools::halley_iterate(f_plus_args,
       guess, min, max, digits, max_iter);
  } catch (const std::exception& e) {
    throw e;
  }
  return ret;
}
*/

// Non-tol versions
template <typename F, typename FDiv, typename FDivDiv, typename GuessScalar, typename MinScalar, typename MaxScalar, typename... Types>
auto root_finder(F&& f, FDiv&& f_div, FDivDiv&& f_div_div, GuessScalar guess, MinScalar min, MaxScalar max, Types&&... args) {
    int digits = 16;
    std::uintmax_t max_iter = 100;
    return root_finder_tol(f, f_div, f_div_div, guess, min, max, digits, max_iter, args...);
}
/*
template <typename F, typename FDiv, typename GuessScalar, typename MinScalar, typename MaxScalar, typename... Types>
double root_finder(F&& f, FDiv&& f_div, GuessScalar guess, MinScalar min, MaxScalar max, Types&&... args) {
    int digits = 16;
    std::uintmax_t max_iter = 100;
    return root_finder_tol(f, f_div, guess, min, max, digits, max_iter, args...);
}

template <typename F, typename GuessScalar, typename MinScalar, typename MaxScalar, typename... Types>
double root_finder(F&& f, GuessScalar guess, MinScalar min, MaxScalar max, Types&&... args) {
    int digits = 16;
    std::uintmax_t max_iter = 100;
    return root_finder_tol(f, guess, min, max, digits, max_iter, args...);
}
*/
}
}
#endif
