#ifndef STAN_MATH_PRIM_FUN_HOLDER_HPP
#define STAN_MATH_PRIM_FUN_HOLDER_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_plain_type.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <memory>
#include <type_traits>
#include <utility>

/**
 * \defgroup eigen_expressions Eigen expressions
 */
/**
 * \ingroup eigen_expressions
 * \defgroup returning_expressions Returning expressions
 */
/**
 * \ingroup returning_expressions
 *
 * Operations in Eigen Expressions hold their arguments either by value or by
 * reference. Which one is chosen depends on type of the argument. Other
 * operations are held by value. "Heavy" objects that can hold data themselves,
 * such as `Eigen::Matrix` or `Eigen::Ref` are instead held by reference. THis
 * is the only criterion - holding rvalue arguments by value is not supported,
 * so we can not use perfect forwarding.
 *
 * When returning an expression from function we have to be careful that any
 * arguments in this expression that are held by reference do not go out of
 * scope. For instance, consider the function:
 *
 * ```
 * template<typename T>
 * auto f(const T& x){
 *   const Eigen::Ref<const Eigen::VectorXd>& x_ref = x;
 *   return x_ref.cwiseProduct(x_ref);
 * }
 * ```
 * And the code calling it:
 * ```
 * Eigen::MatrixXd test_mat(2,2);
 * test_mat << 5, 5, 5, 5;
 * VectorXd X  = f(test_mat.diagonal());
 * ```
 * This function will return back a `CwiseBinaryOp` Eigen expression, which is
 * then evaluated out of the function scope when assigned to `X`. The expression
 * references `x_ref`, which was created withing function and destroyed when
 * the function returned. The returned expression is evaluated after the
 * function returned, so its evaluation references a matrix that was already
 * deleted. In other words the returned expression contains a dangling
 * reference.
 *
 * So a function returning an expression referencing local matrices or
 * matrices that were rvalue reference arguments to the function will not work.
 *
 * A workarount to this issue is allocating and constructing or moving such
 * objects to heap. `Holder` object is a no-op operation that can also take
 * pointers to such objects and release them when it goes out of scope. It can
 * be created either by directly supplying pointers to such objects to `holder`
 * function or by forwarding function arguments and moving local variables to
 * `make_holder`, which will move any rvalues to heap first.
 */

// This was implenmented following the tutorial on edding new expressions to
// Eigen: https://eigen.tuxfamily.org/dox/TopicNewExpressionType.html

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
    // Possible flags are documented here:
    // https://eigen.tuxfamily.org/dox/group__flags.html
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

  // we need to explicitely default copy and move constructors as we are
  // defining copy and move assignment operators
  Holder(const Holder<ArgType, Ptrs...>&) = default;
  Holder(Holder<ArgType, Ptrs...>&&) = default;

  // all these functions just call the same on the argument
  Eigen::Index rows() const { return m_arg.rows(); }
  Eigen::Index cols() const { return m_arg.cols(); }
  Eigen::Index innerStride() const { return m_arg.innerStride(); }
  Eigen::Index outerStride() const { return m_arg.outerStride(); }
  auto* data() { return m_arg.data(); }

  /**
   * Assignment operator assigns expresssions.
   * @param other expression to assign  to this
   * @return *this
   */
  template <typename T, require_eigen_t<T>* = nullptr>
  inline Holder<ArgType, Ptrs...>& operator=(const T& other) {
    m_arg = other;
    return *this;
  }

  // copy and move assignment operators need to be separately overloaded,
  // otherwise defaults will be used.
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
  typedef typename remove_all<ArgType>::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  typedef typename XprType::Scalar Scalar;
  enum {
    CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
    // Possible flags are documented here:
    // https://eigen.tuxfamily.org/dox/group__flags.html
    Flags = evaluator<ArgTypeNestedCleaned>::Flags,
    Alignment = evaluator<ArgTypeNestedCleaned>::Alignment,
  };

  evaluator<ArgTypeNestedCleaned> m_argImpl;

  explicit evaluator(const XprType& xpr) : m_argImpl(xpr.m_arg) {}

  // all these functions just call the same on the argument
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
    return m_argImpl.template writePacket<StoreMode, PacketType>(row, col, x);
  }
  template <int StoreMode, typename PacketType>
  EIGEN_STRONG_INLINE void writePacket(Index index, const PacketType& x) {
    return m_argImpl.template writePacket<StoreMode, PacketType>(index, x);
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
 * Handles single element (moving rvalue non-expressions to heap) for
 * construction of `holder` or `holder_cl` from a functor. For lvalues or
 * expressions just sets the `res` pointer.
 * @tparam T type of the element
 * @param a element to handle
 * @param res resulting pointer to element
 * @return tuple of pointers allocated on heap (empty).
 */
template <typename T>
auto holder_handle_element(T& a, T*& res) {
  res = &a;
  return std::make_tuple();
}
template <typename T,
          std::enable_if_t<!(Eigen::internal::traits<std::decay_t<T>>::Flags
                             & Eigen::NestByRefBit)>* = nullptr>
auto holder_handle_element(T&& a, std::remove_reference_t<T>*& res) {
  res = &a;
  return std::make_tuple();
}

/**
 * Handles single element (moving rvalue non-expressions to heap) for
 * construction of `holder` or `holder_cl` from a functor. Rvalue non-expression
 * is moved to heap and the pointer to heap memory is assigned to res and
 * returned in a tuple.
 * @tparam T type of the element
 * @param a element to handle
 * @param res resulting pointer to element
 * @return tuple of pointers allocated on heap (containing single pointer).
 */
template <typename T, require_t<std::is_rvalue_reference<T&&>>* = nullptr,
          std::enable_if_t<
              static_cast<bool>(Eigen::internal::traits<std::decay_t<T>>::Flags&
                                    Eigen::NestByRefBit)>* = nullptr>
auto holder_handle_element(T&& a, T*& res) {
  res = new T(std::move(a));
  return std::make_tuple(res);
}
template <typename T, require_t<std::is_rvalue_reference<T&&>>* = nullptr,
          require_not_eigen_t<T>* = nullptr>
auto holder_handle_element(T&& a, T*& res) {
  res = new T(std::move(a));
  return std::make_tuple(res);
}

/**
 * Second step in implementation of construction `holder` from a functor.
 * Constructs holder object form given expression and tuple of pointers.
 * @tparam T type of the result expression
 * @tparam Is index sequence for `ptrs`
 * @tparam Args types of pointes to heap
 * @param expr result expression
 * @param ptrs pointers to heap that need to be released when the expression is
 * destructed
 * @return `holder` referencing given expression
 */
template <typename T, std::size_t... Is, typename... Args>
auto make_holder_impl_construct_object(T&& expr, std::index_sequence<Is...>,
                                       const std::tuple<Args*...>& ptrs) {
  return holder(std::forward<T>(expr), std::get<Is>(ptrs)...);
}

/**
 * Implementation of construction `holder` from a functor.
 * @tparam F type of the functor
 * @tparam Is index sequence for `args`
 * @tparam Args types of the arguments
 * @param func functor
 * @param args arguments for the functor
 * @return `holder` referencing expression constructed by given functor
 */
template <typename F, std::size_t... Is, typename... Args>
auto make_holder_impl(const F& func, std::index_sequence<Is...>,
                      Args&&... args) {
  std::tuple<std::remove_reference_t<Args>*...> res;
  auto ptrs = std::tuple_cat(
      holder_handle_element(std::forward<Args>(args), std::get<Is>(res))...);
  return make_holder_impl_construct_object(
      func(*std::get<Is>(res)...),
      std::make_index_sequence<std::tuple_size<decltype(ptrs)>::value>(), ptrs);
}

}  // namespace internal

/**
 * Constructs an expression from given arguments using given functor.
 * This is similar to calling the functor with given arguments. Except that any
 * rvalue argument will be moved to heap first. The arguments moved to heap are
 * deleted once the expression is destructed.
 *
 * @tparam F type of the functor
 * @tparam Args types of the arguments
 * @param func the functor
 * @param args arguments for the functor
 * @return `holder` referencing expression constructed by given functor
 */
template <typename F, typename... Args,
          require_not_plain_type_t<
              decltype(std::declval<F>()(std::declval<Args&>()...))>* = nullptr>
auto make_holder(const F& func, Args&&... args) {
  return internal::make_holder_impl(func,
                                    std::make_index_sequence<sizeof...(Args)>(),
                                    std::forward<Args>(args)...);
}

/**
 * Calls given function with given arguments. No `holder` is necessary if the
 * function is not returning Eigen expression.
 *
 * @tparam F type of the functor
 * @tparam Args types of the arguments
 * @param func the functor
 * @param args arguments for the functor
 * @return `holder` referencing expression constructed by given functor
 */
template <typename F, typename... Args,
          require_plain_type_t<
              decltype(std::declval<F>()(std::declval<Args&>()...))>* = nullptr>
auto make_holder(const F& func, Args&&... args) {
  return func(std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
