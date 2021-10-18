#ifndef STAN_MATH_PRIM_PLUGINS_ADJ_VIEW_H
#define STAN_MATH_PRIM_PLUGINS_ADJ_VIEW_H

template <typename Scalar, typename Enable = void>
struct adj_impl { };

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

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const val_return_t<Scalar>& adj(const Scalar& x) {
  return adj_impl<Scalar>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline val_return_t<Scalar>& adj_ref(Scalar& x) {
  return adj_impl<Scalar>::run(x);
}

template <typename Scalar>
struct scalar_adj_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_adj_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const val_return_t<Scalar>& operator() (const Scalar& a) const {
    return adj(a);
  }
};

template <typename Scalar>
struct scalar_adj_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_adj_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE val_return_t<Scalar>& operator() (const Scalar& a) const {
    return adj_ref(*const_cast<Scalar*>(&a));
  }
};

template <typename T, typename Enable = void>
struct adj_stride { };

template <typename T>
struct adj_stride<T, std::enable_if_t<!is_vari<T>::value
                                      && !is_var<T>::value>> {
  static constexpr int stride = -1;
};

template <typename T>
struct adj_stride<T, std::enable_if_t<is_vari<T>::value>> {
  using vari_t = std::remove_pointer_t<T>;
  static constexpr int stride = sizeof(vari_t) / sizeof(typename vari_t::value_type);
};

template <typename T>
struct adj_stride<T, std::enable_if_t<is_var<T>::value>> {
  using vari_t = typename std::decay_t<std::remove_pointer_t<T>>::vari_type;
  static constexpr int stride = sizeof(vari_t) / sizeof(typename vari_t::value_type);
};

typedef CwiseUnaryOp<scalar_adj_op<Scalar>, const Derived> AdjReturnType;
typedef CwiseUnaryView<scalar_adj_ref_op<Scalar>, Derived,
                       adj_stride<Scalar>::stride,
                       adj_stride<Scalar>::stride> NonConstAdjReturnType;

EIGEN_DEVICE_FUNC
inline const AdjReturnType
adj() const { return AdjReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstAdjReturnType
adj() { return NonConstAdjReturnType(derived()); }

#endif
