#ifndef STAN_MATH_PRIM_PLUGINS_VAL_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_VAL_VIEW_H

template <typename Scalar, typename Enable = void>
struct val_impl { };

template <typename Scalar>
struct val_impl<Scalar, std::enable_if_t<std::is_arithmetic<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x;
  }
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x;
  }
};

template <typename Scalar>
struct val_impl<Scalar, std::enable_if_t<is_vari<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x->val_;
  }
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x->val_;
  }
};

template <typename Scalar>
struct val_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x.vi_->val_;
  }
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x.vi_->val_;
  }
};

template <typename Scalar>
struct val_impl<Scalar, std::enable_if_t<is_fvar<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return x.val_;
  }
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return x.val_;
  }
};

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const val_return_t<Scalar>& val(const Scalar& x) {
  return val_impl<Scalar>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline
val_return_t<Scalar>&
val_ref(Scalar& x) {
  return val_impl<Scalar>::run(x);
}

template<typename Scalar>
struct scalar_val_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_val_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  const val_return_t<Scalar>& operator() (const Scalar& a) const {
    return val(a);
  }
};

template<typename Scalar>
struct scalar_val_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_val_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE   
  val_return_t<Scalar>&
  operator() (const Scalar& a) const {
    return val_ref(*const_cast<Scalar*>(&a));
  }
};

typedef CwiseUnaryOp<scalar_val_op<Scalar>, const Derived> valReturnType;
typedef std::conditional_t<is_var<Scalar>::value || is_vari<Scalar>::value,
                   const valReturnType,
                   CwiseUnaryView<scalar_val_ref_op<Scalar>, Derived>>
NonConstvalReturnType;

EIGEN_DEVICE_FUNC
inline const valReturnType
val() const { return valReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstvalReturnType
val() { return NonConstvalReturnType(derived()); }

#endif
