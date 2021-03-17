#ifndef STAN_MATH_PRIM_FUN_HOLDER_HPP
#define STAN_MATH_PRIM_FUN_HOLDER_HPP

#include <stan/math/prim/meta/index_apply.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_plain_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
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
 * Operations in Eigen Expressions hold their Args either by value or by
 * reference. Which one is chosen depends on type of the argument. "Heavy"
 * objects that can hold data themselves, such as `Eigen::Matrix` or
 * `Eigen::Ref` are held by reference. Other operations are held by value. This
 * is the only criterion - holding rvalue Args by value is not supported, so we
 * can not use perfect forwarding.
 *
 * When returning an expression from function we have to be careful that any
 * Args in this expression that are held by reference do not go out of
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
 * matrices that were rvalue reference Args to the function will not work.
 *
 * A workarount to this issue is any referenced object part of the returned
 * expression. `Holder` object is a no-op operation that can also contain such
 * objects. It can be created either by forwarding function argumnets and moving
 * local variables to `make_holder`, which will construct a Holder object
 * containing these objects and behaving as an expression returned by given
 * functor.
 *
 * Same holder object can also be used for OpenCL kernel generator expressions.
 * Kernel generator expressions can themselves hold rvalue argumets, so using
 * holder for simple cases, such as the example above is not needed. However, we
 * do need holder when a rvalue object is referenced multiple times by the
 * returned expression.
 */

// This was implenmented following the tutorial on edding new expressions to
// Eigen: https://eigen.tuxfamily.org/dox/TopicNewExpressionType.html

namespace stan {
namespace math {

template <class BaseExpr, typename... Ptrs>
class Holder;

}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

template <class BaseExpr, typename... Args>
struct traits<stan::math::Holder<BaseExpr, Args...>> {
  typedef typename BaseExpr::StorageKind StorageKind;
  typedef typename traits<BaseExpr>::XprKind XprKind;
  typedef typename BaseExpr::StorageIndex StorageIndex;
  typedef typename BaseExpr::Scalar Scalar;
  enum {
    // Possible flags are documented here:
    // https://eigen.tuxfamily.org/dox/group__flags.html
    Flags = (BaseExpr::Flags
             & (RowMajorBit | LvalueBit | LinearAccessBit | DirectAccessBit
                | PacketAccessBit | NoPreferredStorageOrderBit))
            | NestByRefBit,
    RowsAtCompileTime = BaseExpr::RowsAtCompileTime,
    ColsAtCompileTime = BaseExpr::ColsAtCompileTime,
    MaxRowsAtCompileTime = BaseExpr::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = BaseExpr::MaxColsAtCompileTime,
    InnerStrideAtCompileTime = BaseExpr::InnerStrideAtCompileTime,
    OuterStrideAtCompileTime = BaseExpr::OuterStrideAtCompileTime
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
template <typename BaseExpr, typename... Args>
class Holder
    : private std::tuple<Args...>,
      public Eigen::internal::dense_xpr_base<Holder<BaseExpr, Args...>>::type {
 public:
  typedef
      typename Eigen::internal::ref_selector<Holder<BaseExpr, Args...>>::type
          Nested;
  typename Eigen::internal::ref_selector<BaseExpr>::non_const_type m_arg;

  template <typename Functor, require_same_t<decltype(std::declval<Functor>()(
                                                 std::declval<Args&>()...)),
                                             BaseExpr>* = nullptr>
  explicit Holder(const Functor& f, Args&&... args)
      : std::tuple<Args...>(std::forward<Args>(args)...),
        m_arg(index_apply<sizeof...(Args)>([this, &f](auto... Is) {
          return f(std::get<Is>(*static_cast<std::tuple<Args...>*>(this))...);
        })) {}

  // we need to explicitely default copy and move constructors as we are
  // defining copy and move assignment operators
  Holder(const Holder<BaseExpr, Args...>&) = default;
  Holder(Holder<BaseExpr, Args...>&&) = default;

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
  inline Holder<BaseExpr, Args...>& operator=(const T& other) {
    m_arg = other;
    return *this;
  }

  // copy and move assignment operators need to be separately overloaded,
  // otherwise defaults will be used.
  inline Holder<BaseExpr, Args...>& operator=(
      const Holder<BaseExpr, Args...>& other) {
    m_arg = other;
    return *this;
  }
  inline Holder<BaseExpr, Args...>& operator=(
      Holder<BaseExpr, Args...>&& other) {
    m_arg = std::move(other.m_arg);
    return *this;
  }
};

}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

template <typename BaseExpr, typename... Args>
struct evaluator<stan::math::Holder<BaseExpr, Args...>>
    : evaluator_base<stan::math::Holder<BaseExpr, Args...>> {
  typedef stan::math::Holder<BaseExpr, Args...> XprType;
  typedef typename remove_all<BaseExpr>::type ArgTypeNestedCleaned;
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
auto make_holder(const F& func, Args&&... args) {
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
auto make_holder(const F& func, Args&&... args) {
  return Holder<std::decay_t<decltype(func(args...))>, Args...>(
      func, std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
