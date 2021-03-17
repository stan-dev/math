#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_holder_cl__HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_holder_cl__HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <memory>
#include <string>
#include <type_traits>
#include <set>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents a no-op in kernel generator expressions. This object also owns
 * pointers to dynamically allocated objects used in its argument expression.
 * When this object is destructed, those objects are deleted.
 * @tparam Derived derived type
 * @tparam T type of the argument
 */
template <typename T, typename... Args>
class holder_cl_
    : private std::tuple<Args...>,
      public operation_cl<holder_cl_<T, Args...>,
                          typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<holder_cl_<T, Args...>, Scalar, T>;
  using base::var_name_;

  /**
   * Constructor
   * @param f functor callable with given arguments that returns an expression
   * @param args arguments for the functor
   */
  template <typename Functor, require_same_t<decltype(std::declval<Functor>()(
                                                 std::declval<Args&>()...)),
                                             T>* = nullptr>
  explicit holder_cl_(const Functor& f, Args&&... args)
      : std::tuple<Args...>(std::forward<Args>(args)...),
        base(index_apply<sizeof...(Args)>([this, &f](auto... Is) {
          return f(std::get<Is>(*static_cast<std::tuple<Args...>*>(this))...);
        })) {}
};

/**
 * Calls given function with given Args. No `holder` is necessary if the
 * function is not returning Eigen expression or if all arguments are lvalues.
 *
 * @tparam F type of the functor
 * @tparam Args types of the Args
 * @param func the functor
 * @param args Args for the functor
 * @return `holder` referencing expression constructed by given functor
 */
template <
    typename F, typename... Args,
    require_t<disjunction<conjunction<std::is_lvalue_reference<Args&&>...>,
                          is_plain_type<decltype(std::declval<F>()(
                              std::declval<Args&>()...))>>>* = nullptr>
auto make_holder_cl(const F& func, Args&&... args) {
  return func(std::forward<Args>(args)...);
}

/**
 * Constructs an expression from given Args using given functor.
 * This is similar to calling the functor with given Args. Except that any
 * rvalue argument will be moved to heap first. The Args moved to heap are
 * deleted once the expression is destructed.
 *
 * @tparam F type of the functor
 * @tparam Args types of the Args
 * @param func the functor
 * @param args Args for the functor
 * @return `holder` referencing expression constructed by given functor
 */
template <typename F, typename... Args,
          require_not_plain_type_t<
              decltype(std::declval<F>()(std::declval<Args&>()...))>* = nullptr,
          require_any_t<std::is_rvalue_reference<Args&&>...>* = nullptr>
auto make_holder_cl(const F& func, Args&&... args) {
  return holder_cl_<std::decay_t<decltype(func(args...))>, Args...>(
      func, std::forward<Args>(args)...);
}

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
