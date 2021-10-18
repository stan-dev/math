#ifndef STAN_MATH_PRIM_PLUGINS_VI_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_VI_VIEW_H

template <typename Scalar, typename Enable = void>
struct vi_impl { };

template <typename Scalar>
struct vi_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline const vi_return_t<Scalar>& run(const Scalar& x) {
    return x.vi_;
  }
  EIGEN_DEVICE_FUNC
  static inline vi_return_t<Scalar>& run(Scalar& x) {
    return x.vi_;
  }
};

template <typename Scalar>
struct scalar_vi_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_vi_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const vi_return_t<Scalar>& operator() (const Scalar& a) const {
    return vi(a);
  }
};

template <typename Scalar>
struct scalar_vi_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_vi_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE vi_return_t<Scalar>& operator() (const Scalar& a) const {
    return vi_ref(*const_cast<Scalar*>(&a));
  }
};

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const vi_return_t<Scalar>& vi(const Scalar& x) {
  return vi_impl<Scalar>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline vi_return_t<Scalar>& vi_ref(Scalar& x) {
  return vi_impl<Scalar>::run(x);
}

typedef CwiseUnaryOp<scalar_vi_op<Scalar>, const Derived> viReturnType;
typedef CwiseUnaryView<scalar_vi_ref_op<Scalar>, Derived> NonConstviReturnType;

EIGEN_DEVICE_FUNC
inline const viReturnType
vi() const { return viReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstviReturnType
vi() { return NonConstviReturnType(derived()); }

#endif
