#ifndef STAN_MATH_PRIM_PLUGINS_D_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_D_VIEW_H


// Struct for returning the gradient from an fvar<T>
template <typename Scalar>
struct d_impl {
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x.d_;
  }
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x.d_;
  }
};

// Struct implementing operator() to be called by CWiseUnaryOp to
//   return the gradient from a const fvar<T>
template <typename Scalar>
struct scalar_d_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_d_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  const val_return_t<Scalar>& operator() (const Scalar& a) const {
    return d_impl<Scalar>::run(a);
  }
};

// Struct implementing operator() to be called by CWiseUnaryView to
//   return the gradient from a non-const fvar<T>
template <typename Scalar>
struct scalar_d_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_d_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  val_return_t<Scalar>& operator() (const Scalar& a) const {
    return d_impl<Scalar>::run(*const_cast<Scalar*>(&a));
  }
};

/**
 * Coefficient-wise function returning a view of the gradients that cannot
 * be modified
 */
EIGEN_DEVICE_FUNC
inline auto d() const {
  return CwiseUnaryOp<scalar_d_op<Scalar>, const Derived>(derived());
}

/**
 * Coefficient-wise function returning a view of the gradients that can
 * be modified.
 */
EIGEN_DEVICE_FUNC
inline auto d() {
  return CwiseUnaryView<scalar_d_ref_op<Scalar>, Derived>(derived());
}

#endif
