#ifndef STAN_MATH_PRIM_FUN_HOLDER_HPP
#define STAN_MATH_PRIM_FUN_HOLDER_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <memory>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

template <class ArgType, typename... Ptrs>
class Holder;

}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

template <class ArgType, typename... Ptrs>
struct traits<stan::math::Holder<ArgType, Ptrs...>> {
  typedef typename ArgType::StorageKind StorageKind;
  typedef typename traits<ArgType>::XprKind XprKind;
  typedef typename ArgType::StorageIndex StorageIndex;
  typedef typename ArgType::Scalar Scalar;
  enum {
    Flags = (ArgType::Flags
             & (RowMajorBit | LvalueBit | LinearAccessBit | DirectAccessBit
                | PacketAccessBit | NoPreferredStorageOrderBit))
            | NestByRefBit,
    RowsAtCompileTime = ArgType::RowsAtCompileTime,
    ColsAtCompileTime = ArgType::ColsAtCompileTime,
    MaxRowsAtCompileTime = ArgType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime,
    InnerStrideAtCompileTime = ArgType::InnerStrideAtCompileTime,
    OuterStrideAtCompileTime = ArgType::OuterStrideAtCompileTime
  };
};

}  // namespace internal
}  // namespace Eigen

namespace stan {
namespace math {

/**
 * A no-op Eigen operation. This object also owns pointers to dynamically
 * allocated objects used in its argument expression. When this object is
 * destructed, those objects are deleted.
 * @tparam Derived derived type
 * @tparam T type of the argument
 */
template <typename ArgType, typename... Ptrs>
class Holder
    : public Eigen::internal::dense_xpr_base<Holder<ArgType, Ptrs...>>::type {
 public:
  typedef typename Eigen::internal::ref_selector<Holder<ArgType, Ptrs...>>::type
      Nested;
  typename Eigen::internal::ref_selector<ArgType>::non_const_type m_arg;
  std::tuple<std::unique_ptr<Ptrs>...> m_unique_ptrs;

  explicit Holder(ArgType&& arg, Ptrs*... pointers)
      : m_arg(arg), m_unique_ptrs(std::unique_ptr<Ptrs>(pointers)...) {}

  Holder(const Holder<ArgType, Ptrs...>&) = default;
  Holder(Holder<ArgType, Ptrs...>&&) = default;

  Eigen::Index rows() const { return m_arg.rows(); }
  Eigen::Index cols() const { return m_arg.cols(); }
  Eigen::Index innerStride() const { return m_arg.innerStride(); }
  Eigen::Index outerStride() const { return m_arg.outerStride(); }

  auto* data() { return m_arg.data(); }

  template <typename T, require_eigen_t<T>* = nullptr>
  inline Holder<ArgType, Ptrs...>& operator=(const T& other) {
    m_arg = other;
    return *this;
  }

  inline Holder<ArgType, Ptrs...>& operator=(
      const Holder<ArgType, Ptrs...>& other) {
    m_arg = other;
    return *this;
  }

  inline Holder<ArgType, Ptrs...>& operator=(Holder<ArgType, Ptrs...>&& other) {
    m_arg = std::move(other.m_arg);
    return *this;
  }
};

}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

template <typename ArgType, typename... Ptrs>
struct evaluator<stan::math::Holder<ArgType, Ptrs...>>
    : evaluator_base<stan::math::Holder<ArgType, Ptrs...>> {
  typedef stan::math::Holder<ArgType, Ptrs...> XprType;
  typedef typename remove_all<typename nested_eval<ArgType, 1>::type>::type
      ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  typedef typename XprType::Scalar Scalar;
  enum {
    CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
    Flags
    = evaluator<ArgTypeNestedCleaned>::Flags
      & (HereditaryBits | LinearAccessBit | PacketAccessBit | DirectAccessBit),
    Alignment = evaluator<ArgTypeNestedCleaned>::Alignment,
  };

  evaluator<ArgTypeNestedCleaned> m_argImpl;

  explicit evaluator(const XprType& xpr) : m_argImpl(xpr.m_arg) {}

  EIGEN_STRONG_INLINE CoeffReturnType coeff(Index row, Index col) const {
    return m_argImpl.coeff(row, col);
  }
  EIGEN_STRONG_INLINE CoeffReturnType coeff(Index index) const {
    return m_argImpl.coeff(index);
  }

  EIGEN_STRONG_INLINE Scalar& coeffRef(Index row, Index col) {
    return m_argImpl.coeffRef(row, col);
  }
  EIGEN_STRONG_INLINE Scalar& coeffRef(Index index) {
    return m_argImpl.coeffRef(index);
  }

  template <int LoadMode, typename PacketType>
  EIGEN_STRONG_INLINE PacketType packet(Index row, Index col) const {
    return m_argImpl.template packet<LoadMode, PacketType>(row, col);
  }
  template <int LoadMode, typename PacketType>
  EIGEN_STRONG_INLINE PacketType packet(Index index) const {
    return m_argImpl.template packet<LoadMode, PacketType>(index);
  }

  template <int StoreMode, typename PacketType>
  EIGEN_STRONG_INLINE void writePacket(Index row, Index col,
                                       const PacketType& x) {
    return m_argImpl.template writePacket(row, col, x);
  }
  template <int StoreMode, typename PacketType>
  EIGEN_STRONG_INLINE void writePacket(Index index, const PacketType& x) {
    return m_argImpl.template writePacket(index, x);
  }
};

}  // namespace internal
}  // namespace Eigen

namespace stan {
namespace math {

/**
 * Constructs a no-op operation that also holds pointer to some other
 * expressions, allocated on heap. When the object is destructed those
 * expressions will be deleted.
 * @tparam T type of argument expression
 * @tparam Ptrs types of pointers
 * @param a argument expression
 * @param ptrs pointers to objects the constructed object will own.
 * @return holder operation
 */
template <typename T, typename... Ptrs,
          std::enable_if_t<sizeof...(Ptrs) >= 1>* = nullptr>
Holder<T, Ptrs...> holder(T&& arg, Ptrs*... pointers) {
  return Holder<T, Ptrs...>(std::forward<T>(arg), pointers...);
}
// trivial case with no pointers constructs no holder object
template <typename T>
T holder(T&& arg) {
  return std::forward<T>(arg);
}

namespace internal {
// the function holder_handle_element is also used in holder_cl
/**
 * Handles single element (moving rvalues to heap) for construction of
 * `holder` or `holder_cl` from a functor. For lvalues just sets the `res`
 * pointer.
 * @tparam T type of the element
 * @param a element to handle
 * @param res resulting pointer to element
 * @return tuple of pointer allocated on heap (empty).
 */
template <typename T>
auto holder_handle_element(T& a, T*& res) {
  res = &a;
  return std::make_tuple();
}

/**
 * Handles single element (moving rvalues to heap) for construction of
 * `holder` or `holder_cl` from a functor. Rvalue is moved to heap and the
 * pointer to heap memory is assigned to res and returned in a tuple.
 * @tparam T type of the element
 * @param a element to handle
 * @param res resulting pointer to element
 * @return tuple of pointer allocated on heap (containing single pointer).
 */
template <typename T>
auto holder_handle_element(std::remove_reference_t<T>&& a, T*& res) {
  res = new T(std::move(a));
  return std::make_tuple(res);
}

/**
 * Second step in implementation of construction `holder` from a functor.
 * @tparam T type of the result expression
 * @tparam Is index sequence for `ptrs`
 * @tparam Args types of pointes to heap
 * @param expr result expression
 * @param ptrs pointers to heap that need to be released when the expression is
 * destructed
 * @return `holder` referencing given expression
 */
template <typename T, std::size_t... Is, typename... Args>
auto make_holder_impl_step2(T&& expr, std::index_sequence<Is...>,
                       const std::tuple<Args*...>& ptrs) {
  return holder(std::forward<T>(expr), std::get<Is>(ptrs)...);
}

/**
 * First step in implementation of construction `holder` from a functor.
 * @tparam T type of the functor
 * @tparam Is index sequence for `args`
 * @tparam Args types of arguments
 * @param func functor
 * @param args arguments for the functor
 * @return `holder` referencing given expression
 */
template <typename T, std::size_t... Is, typename... Args>
auto make_holder_impl_step1(const T& func, std::index_sequence<Is...>,
                      Args&&... args) {
  std::tuple<std::remove_reference_t<Args>*...> res;
  auto ptrs = std::tuple_cat(
      holder_handle_element(std::forward<Args>(args), std::get<Is>(res))...);
  return make_holder_impl_step2(
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
          require_eigen_t<
              decltype(std::declval<T>()(std::declval<Args&>()...))>* = nullptr>
auto make_holder(const T& func, Args&&... args) {
  return internal::make_holder_impl_step1(func,
                                    std::make_index_sequence<sizeof...(Args)>(),
                                    std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
