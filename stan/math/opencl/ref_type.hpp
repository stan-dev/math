#ifndef STAN_MATH_OPENCL_REF_TYPE_HPP
#define STAN_MATH_OPENCL_REF_TYPE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/plain_type.hpp>
#include <stan/math/prim/meta/ref_type.hpp>
#include <type_traits>

namespace stan {

namespace internal {

template <typename T>
struct is_trivial_kg_expression : public std::false_type {};

template <>
struct is_trivial_kg_expression<int> : public std::true_type {};
template <>
struct is_trivial_kg_expression<double> : public std::true_type {};
template <typename T>
struct is_trivial_kg_expression<math::matrix_cl<T>> : public std::true_type {};
template <typename T>
struct is_trivial_kg_expression<math::load_<T>> : public std::true_type {};
template <typename T>
struct is_trivial_kg_expression<math::scalar_<T>> : public std::true_type {};
template <typename T>
struct is_trivial_kg_expression<math::constant_<T>> : public std::true_type {};
template <>
struct is_trivial_kg_expression<math::row_index> : public std::true_type {};
template <>
struct is_trivial_kg_expression<math::col_index> : public std::true_type {};
template <typename T>
struct is_trivial_kg_expression<math::calc_if_<false, T>>
    : public std::true_type {};

template <typename T>
struct is_trivial_kg_expression<math::as_column_vector_or_scalar_<T>>
    : public is_trivial_kg_expression<std::decay_t<T>> {};
template <typename T>
struct is_trivial_kg_expression<math::block_<T>>
    : public is_trivial_kg_expression<std::decay_t<T>> {};
template <typename T, bool Colwise, bool Rowwise>
struct is_trivial_kg_expression<math::broadcast_<T, Colwise, Rowwise>>
    : public is_trivial_kg_expression<std::decay_t<T>> {};
template <typename T>
struct is_trivial_kg_expression<math::calc_if_<true, T>>
    : public is_trivial_kg_expression<std::decay_t<T>> {};
template <typename T>
struct is_trivial_kg_expression<math::holder_cl_<T>>
    : public is_trivial_kg_expression<std::decay_t<T>> {};
template <typename T, bool Colwise, bool Rowwise>
struct is_trivial_kg_expression<math::optional_broadcast_<T, Colwise, Rowwise>>
    : public is_trivial_kg_expression<std::decay_t<T>> {};

}  // namespace internal

template <bool Condition, typename T>
struct ref_type_if<Condition, T, require_all_kernel_expressions_t<T>> {
  using T_plain = plain_type_t<T>;
  using T_optionally_ref
      = std::conditional_t<std::is_rvalue_reference<T>::value,
                           std::remove_reference_t<T>, const T&>;
  using type = std::conditional_t<
      internal::is_trivial_kg_expression<std::decay_t<T>>::value || !Condition,
      T_optionally_ref, T_plain>;
};

}  // namespace stan
#endif
#endif
