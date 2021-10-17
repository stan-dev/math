template <typename T, typename Enable = void>
struct val_return { };

template <typename T>
struct val_return<T, std::enable_if_t<!is_fvar<T>::value>> {
  using type = double;
};
template <typename T>
struct val_return<T, std::enable_if_t<is_fvar<T>::value>> {
  using type = decltype(T::d_);
};

template <typename T>
using val_return_t = typename val_return<std::decay_t<T>>::type;

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline val_return_t<Scalar> val(const Scalar& x){
  return val_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const val_return_t<Scalar>& val_ref(const Scalar& x) {
  return val_ref_impl<Scalar>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline
val_return_t<Scalar>&
val_ref(Scalar& x) {
  return val_ref_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar, typename Enable = void>
struct val_default_impl { };

template <typename Scalar>
struct val_default_impl<Scalar, std::enable_if_t<std::is_arithmetic<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double run(const Scalar& x) {
    return x;
  }
};

template <typename Scalar>
struct val_default_impl<Scalar, std::enable_if_t<is_vari<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double run(const Scalar& x) {
    return x->val_;
  }
};

template <typename Scalar>
struct val_default_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double run(const Scalar& x) {
    return x.vi_->val_;
  }
};

template <typename Scalar>
struct val_default_impl<Scalar, std::enable_if_t<is_fvar<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar> run(const Scalar& x) {
    return x.val_;
  }
};

template<typename Scalar> struct val_impl : val_default_impl<Scalar> {};

template<typename Scalar>
struct scalar_val_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_val_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE val_return_t<Scalar> operator() (const Scalar& a) const { return val(a); }
};

template <typename Scalar, typename Enable = void>
struct val_ref_default_impl { };

template <typename Scalar>
struct val_ref_default_impl<Scalar, std::enable_if_t<std::is_arithmetic<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double& run(Scalar& x) {
    return x;
  }
  EIGEN_DEVICE_FUNC
  static inline const double& run(const Scalar& x) {
    return x;
  }
};

template <typename Scalar>
struct val_ref_default_impl<Scalar, std::enable_if_t<is_fvar<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline val_return_t<Scalar>& run(Scalar& x) {
    return *reinterpret_cast<val_return_t<Scalar>*>(&(x.val_));
  }
  EIGEN_DEVICE_FUNC
  static inline const val_return_t<Scalar>& run(const Scalar& x) {
    return *reinterpret_cast<val_return_t<Scalar>*>(&(x.val_));
  }
};

template <typename Scalar>
struct val_ref_impl : val_ref_default_impl<Scalar> {};

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

template <typename T, typename Enable = void>
struct val_stride { };


template <typename T>
struct val_stride<T, std::enable_if_t<!is_vari<T>::value
                                      && !is_var<T>::value
                                      && !is_fvar<T>::value>> {
  static constexpr int stride = -1;
};

template <typename T>
struct val_stride<T, std::enable_if_t<is_vari<T>::value>> {
  using vari_t = std::remove_pointer_t<T>;
  static constexpr int stride = sizeof(vari_t) / sizeof(typename vari_t::value_type);
};

template <typename T>
struct val_stride<T, std::enable_if_t<is_var<T>::value>> {
  using vari_t = std::remove_pointer_t<decltype(std::declval<T>().vi_)>;
  static constexpr int stride = sizeof(vari_t) / sizeof(typename vari_t::value_type);
};

template <typename T>
struct val_stride<T, std::enable_if_t<is_fvar<T>::value>> {
  using fvar_t = T;
  static constexpr int stride = sizeof(fvar_t) / sizeof(typename fvar_t::Scalar);
};

/** \internal the return type of imag() const */
typedef CwiseUnaryOp<scalar_val_op<Scalar>, const Derived> valReturnType;
/** \internal the return type of imag() */
typedef std::conditional_t<is_var<Scalar>::value || is_vari<Scalar>::value,
                   const valReturnType,
                   CwiseUnaryView<scalar_val_ref_op<Scalar>, Derived,
                                  val_stride<Scalar>::stride,
                                  val_stride<Scalar>::stride>>
NonConstvalReturnType;

EIGEN_DEVICE_FUNC
inline const valReturnType
val() const { return valReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstvalReturnType
val() { return NonConstvalReturnType(derived()); }