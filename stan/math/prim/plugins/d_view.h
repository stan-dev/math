template <typename T>
using d_return_t = decltype(T::d_);

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const d_return_t<Scalar>& d(const Scalar& x) {
  return d_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const d_return_t<Scalar>&
d_ref(const Scalar& x) {
  return d_ref_impl<Scalar>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline d_return_t<Scalar>& d_ref(Scalar& x) {
  return d_ref_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar, typename Enable = void>
struct d_default_impl { };

template <typename Scalar>
struct d_default_impl<Scalar, std::enable_if_t<is_fvar<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline d_return_t<Scalar>& run(Scalar& x) {
    return x.d_;
  }
  EIGEN_DEVICE_FUNC
  static inline const d_return_t<Scalar>& run(const Scalar& x) {
    return x.d_;
  }
};

template <typename Scalar>
struct d_impl : d_default_impl<Scalar> {};

template <typename Scalar>
struct scalar_d_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_d_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE d_return_t<Scalar> operator() (const Scalar& a) const { return d(a); }
};


template <typename Scalar, typename Enable = void>
struct d_ref_default_impl { };

template <typename Scalar>
struct d_ref_default_impl<Scalar, std::enable_if_t<is_fvar<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline d_return_t<Scalar>& run(Scalar& x) {
    return *reinterpret_cast<d_return_t<Scalar>*>(&(x.d_));
  }
  EIGEN_DEVICE_FUNC
  static inline const d_return_t<Scalar>& run(const Scalar& x) {
    return *reinterpret_cast<d_return_t<Scalar>*>(&(x.d_));
  }
};

template <typename Scalar>
struct d_ref_impl : d_ref_default_impl<Scalar> {};

template <typename Scalar>
struct scalar_d_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_d_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE d_return_t<Scalar>& operator() (const Scalar& a) const {
    return d_ref(*const_cast<Scalar*>(&a));
  }
};

template <typename T, typename Enable = void>
struct d_stride {
  static constexpr int stride = 1;
};

template <typename T>
struct d_stride<T, std::enable_if_t<is_fvar<typename T::Scalar>::value>> {
  using fvar_t = std::remove_pointer_t<typename T::Scalar>;
  static constexpr int stride = sizeof(fvar_t) / sizeof(typename fvar_t::Scalar);
};

typedef CwiseUnaryOp<scalar_d_op<Scalar>, const Derived> dReturnType;
typedef CwiseUnaryView<scalar_d_ref_op<Scalar>, Derived,
                       d_stride<Derived>::stride,
                       d_stride<Derived>::stride> NonConstdReturnType;

EIGEN_DEVICE_FUNC
inline const dReturnType
d() const { return dReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstdReturnType
d() { return NonConstdReturnType(derived()); }