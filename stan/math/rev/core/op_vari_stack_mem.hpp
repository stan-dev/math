#ifndef STAN_MATH_REV_CORE_OP_VARI_STACK_MEM_HPP
#define STAN_MATH_REV_CORE_OP_VARI_STACK_MEM_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

namespace internal {

// For basic types just T, vari gives vari*, Eigen gives an Eigen::Map
template <typename T, typename = void>
struct op_vari_tuple_type {
  using type = T;
};

template <typename T>
struct op_vari_tuple_type<T, require_vari_t<T>> {
  using type = T;
};

template <typename T>
struct op_vari_tuple_type<T, require_eigen_t<T>> {
  using type = Eigen::Map<typename std::decay_t<T>::PlainObject>;
};

template <typename T>
using op_vari_tuple_t = typename op_vari_tuple_type<T>::type;

// is a type an eigen type with a scalar of double
template <typename T>
struct is_eigen_arith
    : bool_constant<conjunction<std::is_arithmetic<value_type_t<T>>,
                                is_eigen<T>>::value> {};

// Count the number of eigen matrices with arithmetic types
constexpr size_t op_vari_count_dbl(size_t count) { return count; }

template <typename T>
constexpr size_t op_vari_count_dbl(size_t count, T&& x) {
  return count + is_eigen_arith<T>::value;
}

template <typename T, typename... Types>
constexpr size_t op_vari_count_dbl(size_t count, T&& x, Types&&... args) {
  return op_vari_count_dbl(count + is_eigen_arith<T>::value, args...);
}

template <typename T, require_not_eigen_t<T>* = nullptr>
auto make_op_vari(double**& mem, size_t position, T&& x) {
  return std::forward<T>(x);
}

/**
 * Note: This doesn't work!
 * We want position to be equal to the double pointers position, but this uses
 *  the `args` position. So if we had:
 * {var, var, Eigen::Matrix<double>, var, var, Eigen::Matrix<double>}
 * This currently does
 * {0, 1, (2), 3, 4, (5}
 * Where we want
 * {0, 0 , 0, (0), 0, (1)}
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
auto make_op_vari(double**& mem, size_t position, T&& x) {
  mem[position]
      = ChainableStack::instance_->memalloc_.alloc_array<double>(x.size());
  Eigen::Map<typename std::decay_t<T>::PlainObject>(mem[position], x.rows(),
                                                    x.cols())
      = x;
  return Eigen::Map<typename std::decay_t<T>::PlainObject>(mem[position],
                                                           x.rows(), x.cols());
}

template <size_t... I, typename... Types>
auto make_op_vari_tuple_impl(double**& mem, std::index_sequence<I...>,
                             Types&&... args) {
  return std::make_tuple(make_op_vari(mem, I, args)...);
}

template <typename... Types>
auto make_op_vari_tuple(double**& mem, Types&&... args) {
  return make_op_vari_tuple_impl(mem, std::index_sequence_for<Types...>{},
                                 args...);
}
}  // namespace internal
/**
 * Holds the elements needed in vari operations for the reverse pass and chain
 * call.
 *
 * @tparam Types The types of the operation.
 */
template <typename T, typename... Types>
class op_vari : public vari_value<T> {
 protected:
  double** dbl_mem_;  // Holds mem for eigen matrices of doubles
  std::tuple<internal::op_vari_tuple_t<Types>...>
      vi_;  // Holds the objects needed in the reverse pass.

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
   *
   * Here I think the steps are:
   * 1. Count the number of double matrices and allocate double* pointers for
   * them.
   * 2. Apply a tuple constructor that for vari types just returns itself
   *  and for Eigen matrices of doubles allocates the mem for it on our stack,
   *   then constructs and fills the map.
   */
  op_vari(T val, Types... args)
      : vari_value<T>(val),
        dbl_mem_(ChainableStack::instance_->memalloc_.alloc_array<double*>(
            internal::op_vari_count_dbl(0, args...))),
        vi_(internal::make_op_vari_tuple(dbl_mem_, args...)) {}
};

}  // namespace math
}  // namespace stan
#endif
