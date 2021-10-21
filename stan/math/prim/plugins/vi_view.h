#ifndef STAN_MATH_PRIM_PLUGINS_VI_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_VI_VIEW_H

// Struct for returning the vari* from a var
template <typename Scalar>
struct vi_impl{
  EIGEN_DEVICE_FUNC
  static inline const vi_return_t<Scalar>& run(const Scalar& x) {
    return x.vi_;
  }
  EIGEN_DEVICE_FUNC
  static inline vi_return_t<Scalar>& run(Scalar& x) {
    return x.vi_;
  }
};

// Struct implementing operator() to be called by CWiseUnaryOp to
//   return the vari from a const var
template <typename Scalar>
struct scalar_vi_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_vi_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const vi_return_t<Scalar>& operator() (const Scalar& a) const {
    return vi_impl<Scalar>::run(a);
  }
};

// Struct implementing operator() to be called by CWiseUnaryView to
//   return the vari from a non-const var
template <typename Scalar>
struct scalar_vi_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_vi_ref_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  vi_return_t<Scalar>& operator() (const Scalar& a) const {
    return vi_impl<Scalar>::run(*const_cast<Scalar*>(&a));
  }
};

/**
 * Coefficient-wise function returning a view of the values that cannot
 * be modified
 */
EIGEN_DEVICE_FUNC
inline const auto vi() const {
  return CwiseUnaryOp<scalar_vi_op<Scalar>, const Derived>(derived());
}

/**
 * Coefficient-wise function returning a view of the adjoints that can
 * be modified.
 */
EIGEN_DEVICE_FUNC
inline auto vi() {
  return CwiseUnaryView<scalar_vi_ref_op<Scalar>, Derived> (derived());
}

#endif
