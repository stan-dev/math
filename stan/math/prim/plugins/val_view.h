#ifndef STAN_MATH_PRIM_PLUGINS_VAL_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_VAL_VIEW_H

// Forward declaration to allow specialisations
template <typename Scalar, typename Enable = void>
struct val_impl { };

// Struct for returning an arithmetic input argument unchanged
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

// Struct for returning the value from a vari*
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


// Struct for returning the value from a var
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

// Struct for returning the value from an fvar<T>
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

// Struct implementing operator() to be called by CWiseUnaryOp to
//   return the value from a const var or vari*
template <typename Scalar>
struct scalar_val_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_val_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE
  const val_return_t<Scalar>& operator() (const Scalar& a) const {
    return val_impl<Scalar>::run(a);
  }
};

// Struct implementing operator() to be called by CWiseUnaryView to
//   return the value from a non-const var or vari*
template <typename Scalar>
struct scalar_val_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_val_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE   
  val_return_t<Scalar>&
  operator() (const Scalar& a) const {
    return val_impl<Scalar>::run(*const_cast<Scalar*>(&a));
  }
};

typedef CwiseUnaryOp<scalar_val_op<Scalar>, const Derived> valReturnType;

/**
 * Coefficient-wise function returning a view of the values that cannot
 * be modified
 */
EIGEN_DEVICE_FUNC
inline const auto val() const { return valReturnType(derived()); }

/**
 * Coefficient-wise function returning a view of the adjoints that can
 * be modified.
 *
 * The .val() method will always return a const type for var & vari inputs,
 * so CWiseUnaryOp is called instead.
 */
EIGEN_DEVICE_FUNC
inline auto val() {
  return std::conditional_t<is_var<Scalar>::value || is_vari<Scalar>::value,
                  const valReturnType,
                  CwiseUnaryView<scalar_val_ref_op<Scalar>,Derived>>(derived());
}

#endif
