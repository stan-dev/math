#ifndef STAN_MATH_PRIM_FUNCTOR_ROOT_FINDER_HPP
#define STAN_MATH_PRIM_FUNCTOR_ROOT_FINDER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_bounded.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <boost/math/tools/roots.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {
template <bool ReturnDerivs, typename FRootFunc, typename... Args,
          std::enable_if_t<ReturnDerivs>* = nullptr>
inline auto make_root_func(Args&&... args) {
  return [&args...](auto&& x) {
    return std::decay_t<FRootFunc>::template run<ReturnDerivs>(x, args...);
  };
}

template <bool ReturnDerivs, typename FRootFunc,
          std::enable_if_t<!ReturnDerivs>* = nullptr>
inline auto make_root_func() {
  return [](auto&&... args) {
    return std::decay_t<FRootFunc>::template run<ReturnDerivs>(args...);
  };
}

struct NewtonRootSolver {
  template <typename... Types>
  static inline auto run(Types&&... args) {
    return boost::math::tools::newton_raphson_iterate(
        std::forward<Types>(args)...);
  }
};

struct HalleyRootSolver {
  template <typename... Types>
  static inline auto run(Types&&... args) {
    return boost::math::tools::halley_iterate(std::forward<Types>(args)...);
  }
};

struct SchroderRootSolver {
  template <typename... Types>
  static inline auto run(Types&&... args) {
    return boost::math::tools::schroder_iterate(std::forward<Types>(args)...);
  }
};

}  // namespace internal

/**
 * Solve for root using Boost's Halley method
 * @tparam FRootFunc A struct or class with a static function called `run`.
 *  The structs `run` function must have a boolean template parameter that
 *  when `true` returns a tuple containing the function result and the
 * derivatives needed for the root finder. When the boolean template parameter
 * is `false` the function should return a single value containing the function
 * result.
 * @tparam SolverFun One of the three struct types used to call the root solver.
 *  (`NewtonRootSolver`, `HalleyRootSolver`, `SchroderRootSolver`).
 * @tparam GuessScalar Scalar type
 * @tparam MinScalar Scalar type
 * @tparam MaxScalar Scalar type
 * @tparam Types Arg types to pass to functors in `f_tuple`
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
template <typename FRootFunc, typename SolverFun, typename GuessScalar,
          typename MinScalar, typename MaxScalar, typename... Types,
          require_all_not_st_var<GuessScalar, MinScalar, MaxScalar,
                                 Types...>* = nullptr>
auto root_finder_tol(const GuessScalar guess, const MinScalar min,
                     const MaxScalar max, const int digits,
                     std::uintmax_t& max_iter, Types&&... args) {
  check_bounded("root_finder", "initial guess", guess, min, max);
  check_positive("root_finder", "digits", digits);
  check_positive("root_finder", "max_iter", max_iter);
  using ret_t = return_type_t<GuessScalar, MinScalar, MaxScalar, Types...>;
  ret_t ret = 0;
  auto f_plus_div
      = internal::make_root_func<true, FRootFunc>(std::forward<Types>(args)...);
  try {
    ret = std::decay_t<SolverFun>::run(f_plus_div, ret_t(guess), ret_t(min),
                                       ret_t(max), digits, max_iter);
  } catch (const std::exception& e) {
    throw e;
  }
  return ret;
}

template <typename FRootFunc, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types>
auto root_finder_halley_tol(const GuessScalar guess, const MinScalar min,
                            const MaxScalar max, const int digits,
                            std::uintmax_t& max_iter, Types&&... args) {
  return root_finder_tol<FRootFunc, internal::HalleyRootSolver>(
      guess, min, max, digits, max_iter, std::forward<Types>(args)...);
}

template <typename FRootFunc, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types>
auto root_finder_newton_raphson_tol(const GuessScalar guess,
                                    const MinScalar min, const MaxScalar max,
                                    const int digits, std::uintmax_t& max_iter,
                                    Types&&... args) {
  return root_finder_tol<FRootFunc, internal::NewtonRootSolver>(
      guess, min, max, digits, max_iter, std::forward<Types>(args)...);
}

template <typename FRootFunc, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types>
auto root_finder_schroder_tol(const GuessScalar guess, const MinScalar min,
                              const MaxScalar max, const int digits,
                              std::uintmax_t& max_iter, Types&&... args) {
  return root_finder_tol<FRootFunc, internal::SchroderRootSolver>(
      guess, min, max, digits, max_iter, std::forward<Types>(args)...);
}

/**
 * Solve for root with default values for the tolerances
 * @tparam FRootFunc A struct or class with a static function called `run`.
 *  The structs `run` function must have a boolean template parameter that
 *  when `true` returns a tuple containing the function result and the
 * derivatives needed for the root finder. When the boolean template parameter
 * is `false` the function should return a single value containing the function
 * result.
 * @tparam SolverFun One of the three struct types used to call the root solver.
 *  (`NewtonRootSolver`, `HalleyRootSolver`, `SchroderRootSolver`).
 * @tparam GuessScalar Scalar type
 * @tparam MinScalar Scalar type
 * @tparam MaxScalar Scalar type
 * @tparam Types Arg types to pass to functors in `f_tuple`
 * @param guess An initial guess at the root value
 * @param min The minimum possible value for the result, this is used as an
 * initial lower bracket
 * @param max The maximum possible value for the result, this is used as an
 * initial upper bracket
 * @param args Parameter pack of arguments to pass the the functors in `f_tuple`
 */
template <typename FRootFunc, typename SolverFun, typename GuessScalar,
          typename MinScalar, typename MaxScalar, typename... Types>
auto root_finder(const GuessScalar guess, const MinScalar min,
                 const MaxScalar max, Types&&... args) {
  constexpr int digits = 16;
  std::uintmax_t max_iter = std::numeric_limits<std::uintmax_t>::max();
  return root_finder_tol<FRootFunc, SolverFun>(
      guess, min, max, digits, max_iter, std::forward<Types>(args)...);
}

template <typename FRootFunc, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types>
auto root_finder_hailey(const GuessScalar guess, const MinScalar min,
                        const MaxScalar max, Types&&... args) {
  constexpr int digits = 16;
  std::uintmax_t max_iter = std::numeric_limits<std::uintmax_t>::max();
  return root_finder_halley_tol<FRootFunc>(guess, min, max, digits, max_iter,
                                           std::forward<Types>(args)...);
}

template <typename FRootFunc, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types>
auto root_finder_newton_raphson(const GuessScalar guess, const MinScalar min,
                                const MaxScalar max, Types&&... args) {
  constexpr int digits = 16;
  std::uintmax_t max_iter = std::numeric_limits<std::uintmax_t>::max();
  return root_finder_newton_raphson_tol<FRootFunc>(
      guess, min, max, digits, max_iter, std::forward<Types>(args)...);
}

template <typename FRootFunc, typename GuessScalar, typename MinScalar,
          typename MaxScalar, typename... Types>
auto root_finder_schroder(const GuessScalar guess, const MinScalar min,
                          const MaxScalar max, Types&&... args) {
  constexpr int digits = 16;
  std::uintmax_t max_iter = std::numeric_limits<std::uintmax_t>::max();
  return root_finder_schroder_tol<FRootFunc>(guess, min, max, digits, max_iter,
                                             std::forward<Types>(args)...);
}

}  // namespace math
}  // namespace stan
#endif
