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
namespace internal {
template <typename Tuple, typename... Args>
inline auto func_with_derivs(Tuple&& f_tuple, Args&&... args) {
  return stan::math::apply(
      [&args...](auto&&... funcs) {
        return [&args..., &funcs...](auto&& g) {
          return std::make_tuple(funcs(g, args...)...);
        };
      },
      f_tuple);
}
}  // namespace internal

/**
 * Solve for root using Boost's Halley method
 * @tparam FTuple A tuple holding functors whose signatures all match
 *  `(GuessScalar g, Types&&... Args)`.
 * @tparam GuessScalar Scalar type
 * @tparam MinScalar Scalar type
 * @tparam MaxScalar Scalar type
 * @tparam Types Arg types to pass to functors in `f_tuple`
 * @param f_tuple A tuple of functors to calculate the value and any derivates
 * needed.
 * @param guess An initial guess at the root value
 * @param min The minimum possible value for the result, this is used as an
 * initial lower bracket
 * @param max The maximum possible value for the result, this is used as an
 * initial upper bracket
 * @param digits The desired number of binary digits precision
 * @param max_iter An optional maximum number of iterations to perform. On exit,
 * this is updated to the actual number of iterations performed
 * @param args Parameter pack of arguments to pass the the functors in `f_tuple`
 */
template <typename FTuple, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types,
          require_all_not_st_var<GuessScalar, MinScalar, MaxScalar,
                                 Types...>* = nullptr>
auto root_finder_tol(FTuple&& f_tuple, const GuessScalar guess, const MinScalar min,
                     const MaxScalar max, const int digits, std::uintmax_t& max_iter,
                     Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  check_positive("root_finder", "digits", digits);
  check_positive("root_finder", "max_iter", max_iter);
  using ret_t = return_type_t<GuessScalar, MinScalar, MaxScalar, Types...>;
  ret_t ret = 0;
  auto f_plus_div = internal::func_with_derivs(f_tuple, args...);
  try {
    ret = boost::math::tools::halley_iterate(f_plus_div, ret_t(guess), ret_t(min), ret_t(max),
                                             digits, max_iter);
  } catch (const std::exception& e) {
    throw e;
  }
  return ret;
}

/**
 * Solve for root using Boost Halley method with default values for the
 * tolerances
 * @tparam FTuple A tuple holding functors whose signatures all match
 *  `(GuessScalar g, Types&&... Args)`.
 * @tparam GuessScalar Scalar type
 * @tparam MinScalar Scalar type
 * @tparam MaxScalar Scalar type
 * @tparam Types Arg types to pass to functors in `f_tuple`
 * @param f_tuple A tuple of functors to calculate the value and any derivates
 * needed.
 * @param guess An initial guess at the root value
 * @param min The minimum possible value for the result, this is used as an
 * initial lower bracket
 * @param max The maximum possible value for the result, this is used as an
 * initial upper bracket
 * @param args Parameter pack of arguments to pass the the functors in `f_tuple`
 */
template <typename FTuple, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types>
auto root_finder(FTuple&& f_tuple, const GuessScalar guess, const MinScalar min,
                 const MaxScalar max, Types&&... args) {
  constexpr int digits = 16;
  std::uintmax_t max_iter = std::numeric_limits<std::uintmax_t>::max();
  return root_finder_tol(std::forward<FTuple>(f_tuple), guess, min, max, digits,
                         max_iter, std::forward<Types>(args)...);
}
}  // namespace math
}  // namespace stan
#endif
