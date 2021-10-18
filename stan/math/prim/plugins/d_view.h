#ifndef STAN_MATH_PRIM_PLUGINS_D_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_D_VIEW_H

template <typename Scalar, typename Enable = void>
struct d_impl { };

template <typename Scalar>
struct d_impl<Scalar, std::enable_if_t<is_fvar<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x.d_;
  }
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x.d_;
  }
};

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const val_return_t<Scalar>& d(const Scalar& x) {
  return d_impl<Scalar>::run(x);
}
template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline val_return_t<Scalar>& d_ref(Scalar& x) {
  return d_impl<Scalar>::run(x);
}

template <typename Scalar>
struct scalar_d_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_d_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  const val_return_t<Scalar>& operator() (const Scalar& a) const { return d(a); }
};

template <typename Scalar>
struct scalar_d_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_d_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  val_return_t<Scalar>& operator() (const Scalar& a) const {
    return d_ref(*const_cast<Scalar*>(&a));
  }
};

typedef CwiseUnaryOp<scalar_d_op<Scalar>, const Derived> dReturnType;
typedef CwiseUnaryView<scalar_d_ref_op<Scalar>, Derived> NonConstdReturnType;

EIGEN_DEVICE_FUNC
inline const dReturnType
d() const { return dReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstdReturnType
d() { return NonConstdReturnType(derived()); }

#endif
