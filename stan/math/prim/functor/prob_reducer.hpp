#ifndef STAN_MATH_PRIM_FUNCTOR_PROB_REDUCER_HPP
#define STAN_MATH_PRIM_FUNCTOR_PROB_REDUCER_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Used by distributions to decide whether return type shoudl be Scalar or
 * Vector.
 */
enum class ProbReturnType { Scalar, Vector };

/**
 * For scalars performs summations and is a no-op reducer for eigen vectors.
 * Used in the probability distributions for scalar or vector return types.
 */
template <typename T, typename = void>
struct prob_reducer;

/**
 * For scalars performs summations when given eigen types.
 * @tparam A stan scalar type
 */
template <typename T>
struct prob_reducer<T, require_stan_scalar_t<T>> {
  T ret_;  // Underlying return type

  /**
   * Construct from an Eigen type
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   * @param x will be summed and passed into `ret_`.
   */
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  prob_reducer(EigArr&& x) : ret_(sum(std::forward<EigArr>(x))) {}

  /**
   * Construct from an Eigen type while ignoring size argument passed.
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   * @tparam Tossed An integral type
   * @param x will be summed and passed into `ret_`.
   */
  template <typename EigArr, typename Tossed,
            require_eigen_t<EigArr>* = nullptr,
            require_integral_t<Tossed>* = nullptr>
  prob_reducer(EigArr&& x, Tossed&& /* */)
      : ret_(sum(std::forward<EigArr>(x))) {}

  /**
   * Construct from a scalar type.
   * @tparam Scalar a scalar
   * @param x passed to `ret_`.
   */
  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  prob_reducer(Scalar&& x) : ret_(x) {}

  /**
   * Construct from an Eigen type while ignoring size argument passed.
   * @tparam Scalar a scalar
   * @tparam Tossed an integral type
   * @param x will be summed and inserted into `ret_`.
   */
  template <typename Scalar, typename Tossed,
            require_stan_scalar_t<Scalar>* = nullptr>
  prob_reducer(Scalar&& x, Tossed&& /* */) : ret_(x) {}

  /**
   * Perform summation and then assignment
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   */
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator=(EigArr&& x) {
    ret_ = sum(x);
    return *this;
  }

  /**
   * Assignment
   * @tparam Scalar A scalar type
   */
  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator=(Scalar x) {
    ret_ = x;
    return *this;
  }

  /**
   * Perform summation and then `+=`
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   * @param x Eigen object to be summed.
   */
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator+=(EigArr&& x) {
    ret_ += sum(x);
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator+=(Scalar&& x) {
    ret_ += x;
    return *this;
  }

  /**
   * Perform summation and then `-=`
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   * @param x Eigen object to be summed.
   */
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator-=(EigArr&& x) {
    ret_ -= sum(x);
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator-=(Scalar&& x) {
    ret_ -= x;
    return *this;
  }

  /**
   * Return the underlying scalar return type.
   */
  inline auto ret() noexcept { return ret_; }

  /**
   * Return a zero value, used when distribution has special cases that
   *  immedietly return zero.
   * @tparam Types types to deduce the overall return type of the function.
   */
  template <typename... Types>
  static auto zero(int /* */) {
    return return_type_t<Types...>(0);
  }
};

template <typename T>
struct prob_reducer<T, require_eigen_t<T>> {
  T ret_;

  /**
   * Construct from an Eigen type
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   * @param x will be forwarded into `ret_`.
   */
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  prob_reducer(EigArr&& x) : ret_(std::forward<EigArr>(x)) {}

  /**
   * Construct from an Eigen type while ignoring size argument passed.
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   * @tparam Tossed An integral type
   * @param x will be forwarded to `ret_`.
   */
  template <typename EigArr, typename Size, require_eigen_t<EigArr>* = nullptr,
            require_integral_t<Size>* = nullptr>
  prob_reducer(EigArr&& x, Size /* x */) : ret_(std::forward<EigArr>(x)) {}

  /**
   * Construct from a scalar type.
   * @tparam Scalar a scalar
   * @tparam Size An integral type
   * @param x passed to `ret_` along with size to fill with a base value.
   * @param n The size `ret_` should be
   */
  template <typename Scalar, typename Size,
            require_stan_scalar_t<Scalar>* = nullptr,
            require_integral_t<Size>* = nullptr>
  prob_reducer(Scalar constant, Size n) : ret_(T::Constant(n, constant)) {}

  /**
   * assignment
   * @tparam EigArr A type inheriting from `Eigen::EigenBase`
   */
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator=(EigArr&& x) {
    ret_ = std::forward<EigArr>(x);
    return *this;
  }

  /**
   * assignm scalar by propogating value over `ret_`
   * @tparam Scalar a stan scalar
   * @param x The value to fill `ret_` with.
   */
  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator=(Scalar x) {
    ret_ = Eigen::Array<value_type_t<T>, -1, 1>::Constant(x, ret_.size());
    return *this;
  }

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator+=(EigArr&& x) {
    ret_ += std::forward<EigArr>(x);
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator+=(Scalar&& x) {
    ret_ += x;
    return *this;
  }

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator-=(EigArr&& x) {
    ret_ -= std::forward<EigArr>(x);
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator-=(Scalar&& x) {
    ret_ -= x;
    return *this;
  }

  /**
   * Return the underlying scalar return type. Passed the underlying by
   * moving it which can cause `ret_` to be uninitialized after.
   */
  inline auto&& ret() noexcept { return std::move(ret_); }

  /**
   * Return a zero value, used when distribution has special cases that
   *  immedietly return zero.
   * @tparam Types types to deduce the overall return type of the function.
   * @param size The size of the array to return
   */
  template <typename... Types>
  static auto zero(int size) {
    return Eigen::Array<return_type_t<Types...>, -1, 1>::Constant(0, size)
        .eval();
  }
};

/**
 * Generate a reducer with correct return type.
 * @tparam ReturnType Either Scalar or Vector.
 * @tparam Types A parameter pack of types to deduce the underlying scalar type
 * from
 */
template <ProbReturnType ReturnType, typename... Types>
using prob_return_t = prob_reducer<std::conditional_t<
    ReturnType == ProbReturnType::Scalar, return_type_t<Types...>,
    Eigen::Array<return_type_t<Types...>, -1, 1>>>;

}  // namespace math
}  // namespace stan

#endif
