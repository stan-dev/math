#ifndef STAN_MATH_REV_CORE_OP_VARI_HPP
#define STAN_MATH_REV_CORE_OP_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

namespace internal {

// For basic types just T, vari gives vari, Eigen gives an Eigen::Map
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

// is a type an eigen type with an arithmetic scalar
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
decltype(auto) make_op_vari(double**& mem, size_t position, T&& x) {
  return std::forward<T>(x);
}

/**
 * Allocate memory on the stack for Eigen types with arithmetic scalar
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

/**
 * Helper implimentation that parses out the index sequence for figuring out
 * which arguments are eigen types that need memory allocated on the stack.
 * @tparam I a non type template parameter holding used to initialize the index
 *  sequence and expand the array holding the pointer positions for each
 *  eigen arithmetic type in the argument pack.
 * @tparam Types the parameter pack of arguments to construct the tuple from.
 *   `vari_value` and arithmetic types will be passed along while eigen types
 *   with arithmetic scalars have their data copied to the autodiff stack
 *   and are then represented as an `Eigen::Map` whose subtype is the
 *   `PlainObject` type within the Eigen expression.
 * @param mem Pointer to memory from `op_vari` class for arithmetic eigen types
 * @param positions The array holding the positions of allocated memory in the
 *  double pointer of pointers for each eigen type with arithmetic scalars.
 * @param args parameter pack of arguments to initalize for op_vari.
 * Think of the parameter pack expansion below to be similar to
 * ```
 *   forward_as_tuple(make_op_vari(mem, pos[0], args[0]),
 *      make_op_vari(mem, pos[1], args[1]),
 *      make_op_vari(mem, pos[2], args[2]))
 * ```
 */
template <size_t... I, typename... Types>
auto make_op_vari_tuple_impl(
    double**& mem, const std::array<size_t, sizeof...(Types)>& positions,
    std::index_sequence<I...> /* ignore */, Types&&... args) {
  return std::forward_as_tuple(
      make_op_vari(mem, positions[I], std::forward<Types>(args))...);
}

/**
 * Construct the tuple holding the types for the operation.
 * @tparam Types parameter pack of vari_value<T>* pointers eigen types with
 *  arithmetic scalar types.
 * @param mem Pointer to memory from `op_vari` class for arithmetic eigen types
 * @param args parameter pack of arguments to initalize for op_vari.
 */
template <typename... Types>
auto make_op_vari_tuple(double**& mem, Types&&... args) {
  std::array<size_t, sizeof...(Types)> positions_vec{
      {is_eigen_arith<Types>::value...}};
  size_t accum_val = -1;
  for (int i = 0; i < positions_vec.size(); ++i) {
    if (positions_vec[i] != 0) {
      positions_vec[i] += accum_val;
      ++accum_val;
    }
  }
  return make_op_vari_tuple_impl(mem, positions_vec,
                                 std::index_sequence_for<Types...>{}, args...);
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
  const auto& vi() const { return vi_; }
  auto& avi() { return std::get<0>(vi_); }
  const auto& avi() const { return std::get<0>(vi_); }
  auto& ad() { return std::get<0>(vi_); }
  const auto& ad() const { return std::get<0>(vi_); }
  auto& bvi() { return std::get<1>(vi_); }
  const auto& bvi() const { return std::get<1>(vi_); }
  auto& bd() { return std::get<1>(vi_); }
  const auto& bd() const { return std::get<1>(vi_); }
  auto& cvi() { return std::get<2>(vi_); }
  const auto& cvi() const { return std::get<2>(vi_); }
  auto& cd() { return std::get<2>(vi_); }
  const auto& cd() const { return std::get<2>(vi_); }
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
        vi_(internal::make_op_vari_tuple(dbl_mem_,
                                         std::forward<Types>(args)...)) {}
};

}  // namespace math
}  // namespace stan
#endif
