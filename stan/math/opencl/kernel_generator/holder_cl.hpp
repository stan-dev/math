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
template <typename T, typename... Ptrs>
class holder_cl_
    : public operation_cl<holder_cl_<T, Ptrs...>,
                          typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<holder_cl_<T, Ptrs...>, Scalar, T>;
  using base::var_name_;

 protected:
  std::tuple<std::unique_ptr<Ptrs>...> m_unique_ptrs;

 public:
  /**
   * Constructor
   * @param a argument expression
   * @param ptrs pointers this object takes ownership of
   */
  holder_cl_(T&& a, Ptrs*... ptrs)
      : base(std::forward<T>(a)),
        m_unique_ptrs(std::unique_ptr<Ptrs>(ptrs)...) {}
};

/**
 * Constructs a no-op operation that also holds pointer to some other
 * expressions, allocated on heap. When the object is destructed those
 * expressions will be deleted.
 * @tparam T type of argument expression
 * @tparam Ptrs types of pointers
 * @param a argument expression
 * @param ptrs pointers to objects the constructed object will own.
 * @return holder_cl_ operation
 */
template <typename T, typename... Ptrs,
          require_all_kernel_expressions_t<T, Ptrs...>* = nullptr>
auto holder_cl(T&& a, Ptrs*... ptrs) {
  return holder_cl_<as_operation_cl_t<T>, Ptrs...>(
      as_operation_cl(std::forward<T>(a)), ptrs...);
}

namespace internal {
/**
 * Second step in implementation of construction `holder_cl` from a functor.
 * @tparam T type of the result expression
 * @tparam Is index sequence for `ptrs`
 * @tparam Args types of pointes to heap
 * @param expr result expression
 * @param ptrs pointers to heap that need to be released when the expression is
 * destructed
 * @return `holder_cl` referencing given expression
 */
template <typename T, std::size_t... Is, typename... Args>
auto make_holder_cl_impl_step2(T&& expr, std::index_sequence<Is...>,
                               const std::tuple<Args*...>& ptrs) {
  return holder_cl(std::forward<T>(expr), std::get<Is>(ptrs)...);
}

/**
 * First step in implementation of construction `holder_cl` from a functor.
 * @tparam T type of the functor
 * @tparam Is index sequence for `args`
 * @tparam Args types of arguments
 * @param func functor
 * @param args arguments for the functor
 * @return `holder_cl` referencing given expression
 */
template <typename T, std::size_t... Is, typename... Args>
auto make_holder_cl_impl_step1(const T& func, std::index_sequence<Is...>,
                               Args&&... args) {
  std::tuple<std::remove_reference_t<Args>*...> res;
  auto ptrs = std::tuple_cat(
      holder_handle_element(std::forward<Args>(args), std::get<Is>(res))...);
  return make_holder_cl_impl_step2(
      func(*std::get<Is>(res)...),
      std::make_index_sequence<std::tuple_size<decltype(ptrs)>::value>(), ptrs);
}

}  // namespace internal

/**
 * Constructs an expression from given arguments using given functor.
 * This is similar to calling the functor with given arguments. Except that any
 * rvalue argument will be moved to heap first. The arguments moved to heap are
 * deleted once the expression is destructed.
 * @tparam T type of functor
 * @tparam Args types of arguments
 * @param func the functor
 * @param args arguments for the functor
 */
template <typename T, typename... Args,
          require_all_kernel_expressions_t<
              decltype(std::declval<T>()(std::declval<Args&>()...)),
              Args...>* = nullptr>
auto make_holder_cl(const T& func, Args&&... args) {
  return internal::make_holder_cl_impl_step1(
      func, std::make_index_sequence<sizeof...(Args)>(),
      std::forward<Args>(args)...);
}

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
