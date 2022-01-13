#ifndef STAN_MATH_PRIM_PLUGINS_ADJ_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_ADJ_VIEW_H

// Forward declaration to allow specialisations
template <typename Scalar, typename Enable = void>
struct adj_impl { };

// Struct for returning the adjoint from a vari*
template <typename Scalar>
struct adj_impl<Scalar, std::enable_if_t<is_vari<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x->adj_;
  }
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x->adj_;
  }
};

// Struct for returning the adjoint from a var
template <typename Scalar>
struct adj_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x.vi_->adj_;
  }
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x.vi_->adj_;
  }
};

// Struct implementing operator() to be called by CWiseUnaryOp to
//   return the adjoint from a const var or vari*
template <typename Scalar>
struct scalar_adj_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_adj_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const val_return_t<Scalar>& operator() (const Scalar& a) const {
    return adj_impl<Scalar>::run(a);
  }
};

// Struct implementing operator() to be called by CWiseUnaryView to
//   return the adjoint from a non-const var or vari*
template <typename Scalar>
struct scalar_adj_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_adj_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE val_return_t<Scalar>& operator() (const Scalar& a) const {
    return adj_impl<Scalar>::run(*const_cast<Scalar*>(&a));
  }
};

/**
 * Eigen's CWiseUnaryView deduces the Stride between values in memory as
 * the ratio of the sizeof the Matrix scalar type and the returned scalar type.
 * This means that, by default, the stride for var_value<T> is deduced as
 * sizeof(var) / sizeof(T).
 *
 * However, because var types are pointers-to-implementation variables, the stride
 * actually needs to be be sizeof(vari) / sizeof(T). This struct returns the correct
 * stride when var types are used, otherwise it returns the default argument (-1)
 * which indicates that Eigen should calculate the stride as usual
 */
template <typename T, typename Enable = void>
struct adj_stride { };

template <typename T>
struct adj_stride<T, std::enable_if_t<!is_var<T>::value>> {
  static constexpr int stride = -1;
};

template <typename T>
struct adj_stride<T, std::enable_if_t<is_var<T>::value>> {
  using vari_t = typename std::decay_t<std::remove_pointer_t<T>>::vari_type;
  static constexpr int stride = sizeof(vari_t)
                                  / sizeof(typename vari_t::value_type);
};

/**
 * Coefficient-wise function returning a view of the adjoints that cannot
 * be modified
 */
EIGEN_DEVICE_FUNC
inline const auto adj() const {
  return CwiseUnaryOp<scalar_adj_op<Scalar>, const Derived>(derived());
}

/**
 * Coefficient-wise function returning a view of the adjoints that can
 * be modified. The stride is explicitly specified for var types.
 */
EIGEN_DEVICE_FUNC
inline auto adj() {
  return CwiseUnaryView<scalar_adj_ref_op<Scalar>, Derived,
                       adj_stride<Scalar>::stride,
                       adj_stride<Scalar>::stride>(derived());
}

#endif
