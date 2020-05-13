#ifndef STAN_MATH_REV_CORE_OP_VARI_STACK_MEM_HPP
#define STAN_MATH_REV_CORE_OP_VARI_STACK_MEM_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <tuple>

namespace stan {
namespace math {

template <typename T, typename = void>
struct op_vari_tuple_type {
  using type = T;
}

template <typename T>
struct op_vari_tuple_type<T, require_vari_t<T>> {
  using type = T*;
}

template <typename T>
struct op_vari_tuple_type<T, require_eigen_t<T>> {
  using type = Eigen::Map<typename T::PlainObject>;
}

template <typename T>
using op_vari_tuple_t = op_vari_tuple_type<T>::type;

template <typename T, typename = void>
struct op_vari_tuple_mem_type {
  using type = void*;
}

template <typename T>
struct op_vari_tuple_mem_type<T, require_vari_t<T>> {
  using type = void*;
}

template <typename T>
struct op_vari_tuple_mem_type<T, require_eigen_t<T>> {
  using type = value_type_t<T>*;
}

template <typename T>
using op_vari_tuple_mem_t = op_vari_tuple_mem_type<T>::type;

template <typename T>
struct is_eigen_dbl : bool_constant<conjunction<std::is_arithmetic<value_type_t<T>>, is_eigen<T>>::value> {};

constexpr size_t op_vari_count_dbl(size_t count) {
    return count;
}

template <typename T>
constexpr size_t op_vari_count_dbl(size_t count, T&& x) {
    return count + is_eigen_dbl<T>::value;
}

template <typename T, typename... Types>
constexpr size_t op_vari_count_dbl(size_t count, T&& x, Types&&... args) {
    return op_vari_count_dbl(count + is_eigen_dbl<T>::value, args...);
}

template <typename... Types>
auto make_op_vari_tuple(double** mem, Types&&... args) {
  constexpr int num_dbls = std::make_integer_sequence<int, op_vari_count_dbl(0, args...)>;
  return std::make_tuple(make_op_vari(mem, args)...)
}
/**
 * Holds the elements needed in vari operations for the reverse pass and chain
 * call.
 *
 * @tparam Types The types of the operation.
 */
template <typename T, typename... Types>
class op_vari : public vari_value<T> {
 protected:
  double** dbl_mem_; // Holds mem for eigen matrices of doubles
  std::tuple<op_vari_tuple_t<Types>...> vi_;  // Holds the objects needed in the reverse pass.

 public:
  /**
   * Get an element from the tuple of vari ops. Because of name lookup rules
   *  this function needs to be called as \c this->template get<N>()
   * @tparam Ind The index of the tuple to retrieve.
   * @return the element inside of the tuple at index Ind.
   */
  template <std::size_t Ind>
  auto& get() {
    return std::get<Ind>(vi_);
  }

  /**
   * Return a reference to the tuple holding the vari ops. This is commonly
   *  used in conjunction with \c std::get<N>()
   * @return The tuple holding the vari ops.
   */
  auto& vi() { return vi_; }
  auto& avi() { return std::get<0>(vi_); }
  auto& ad() { return std::get<0>(vi_); }
  auto& bvi() { return std::get<1>(vi_); }
  auto& bd() { return std::get<1>(vi_); }
  auto& cvi() { return std::get<2>(vi_); }
  auto& cd() { return std::get<2>(vi_); }
  /**
   * Constructor for passing in vari and ops objects.
   * @param val Value to initialize the vari to.
   * @param args Ops passed into the tuple and used later in chain method.
   */
  op_vari(T val, Types... args)
      : vari_value<T>(val),
      dbl_mem_(ChainableStack::instance_->memalloc_.alloc_array<double*>(op_vari_count_dbl(0, args...))),
      vi_(make_op_vari_tuple(dbl_mem_, args...)) {}
};

}  // namespace math
}  // namespace stan
#endif
