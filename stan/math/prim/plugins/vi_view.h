template <typename T>
using vi_return_t = decltype(T::vi_);

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const vi_return_t<Scalar>& vi(const Scalar& x) {
  return vi_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const vi_return_t<Scalar>&
vi_ref(const Scalar& x) {
  return vi_ref_impl<Scalar>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline vi_return_t<Scalar>& vi_ref(Scalar& x) {
  return vi_ref_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar, typename Enable = void>
struct vi_default_impl { };

template <typename Scalar>
struct vi_default_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline vi_return_t<Scalar>& run(Scalar& x) {
    return x.vi_;
  }
  EIGEN_DEVICE_FUNC
  static inline const vi_return_t<Scalar>& run(const Scalar& x) {
    return x.vi_;
  }
};

template <typename Scalar>
struct vi_impl : vi_default_impl<Scalar> {};

template <typename Scalar>
struct scalar_vi_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_vi_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE vi_return_t<Scalar> operator() (const Scalar& a) const {
    return vi(a);
  }
};


template <typename Scalar, typename Enable = void>
struct vi_ref_default_impl { };

template <typename Scalar>
struct vi_ref_default_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline vi_return_t<Scalar>& run(Scalar& x) {
    return *reinterpret_cast<vi_return_t<Scalar>*>(&(x.vi_));
  }
  EIGEN_DEVICE_FUNC
  static inline const vi_return_t<Scalar>& run(const Scalar& x) {
    return *reinterpret_cast<vi_return_t<Scalar>*>(&(x.vi_));
  }
};

template <typename Scalar>
struct vi_ref_impl : vi_ref_default_impl<Scalar> {};

template <typename Scalar>
struct scalar_vi_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_vi_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE vi_return_t<Scalar>& operator() (const Scalar& a) const {
    return vi_ref(*const_cast<Scalar*>(&a));
  }
};

typedef CwiseUnaryOp<scalar_vi_op<Scalar>, const Derived> viReturnType;
typedef CwiseUnaryView<scalar_vi_ref_op<Scalar>, Derived> NonConstviReturnType;

EIGEN_DEVICE_FUNC
inline const viReturnType
vi() const { return viReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstviReturnType
vi() { return NonConstviReturnType(derived()); }