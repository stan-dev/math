template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const double& adj(const Scalar& x) {
  return adj_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline const double&
adj_ref(const Scalar& x) {
  return adj_ref_impl<Scalar>::run(x);
}

template <typename Scalar>
EIGEN_DEVICE_FUNC
static inline double& adj_ref(Scalar& x) {
  return adj_ref_impl<eigen_base_filter_t<Scalar>>::run(x);
}

template <typename Scalar, typename Enable = void>
struct adj_default_impl { };

template <typename Scalar>
struct adj_default_impl<Scalar, std::enable_if_t<is_vari<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double& run(Scalar& x) {
    return x->adj_;
  }
  EIGEN_DEVICE_FUNC
  static inline const double& run(const Scalar& x) {
    return x->adj_;
  }
};

template <typename Scalar>
struct adj_default_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double& run(Scalar& x) {
    return x.vi_->adj_;
  }
  EIGEN_DEVICE_FUNC
  static inline const double& run(const Scalar& x) {
    return x.vi_->adj_;
  }
};

template <typename Scalar>
struct adj_impl : adj_default_impl<Scalar> {};

template <typename Scalar>
struct scalar_adj_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_adj_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE double operator() (const Scalar& a) const { return adj(a); }
};


template <typename Scalar, typename Enable = void>
struct adj_ref_default_impl { };

template <typename Scalar>
struct adj_ref_default_impl<Scalar, std::enable_if_t<is_vari<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double& run(Scalar& x) {
    return *reinterpret_cast<double*>(&(x->adj_));
  }
  EIGEN_DEVICE_FUNC
  static inline const double& run(const Scalar& x) {
    return *reinterpret_cast<double*>(&(x->adj_));
  }
};

template <typename Scalar>
struct adj_ref_default_impl<Scalar, std::enable_if_t<is_var<Scalar>::value>> {
  EIGEN_DEVICE_FUNC
  static inline double& run(Scalar& x) {
    return *reinterpret_cast<double*>(&(x.vi_->adj_));
  }
  EIGEN_DEVICE_FUNC
  static inline const double& run(const Scalar& x) {
    return *reinterpret_cast<double*>(&(x.vi_->adj_));
  }
};

template <typename Scalar>
struct adj_ref_impl : adj_ref_default_impl<Scalar> {};

template <typename Scalar>
struct scalar_adj_ref_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_adj_ref_op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE double& operator() (const Scalar& a) const {
    return adj_ref(*const_cast<Scalar*>(&a));
  }
};

typedef CwiseUnaryOp<scalar_adj_op<Scalar>, const Derived> AdjReturnType;
typedef CwiseUnaryView<scalar_adj_ref_op<Scalar>, Derived> NonConstAdjReturnType;

EIGEN_DEVICE_FUNC
inline const AdjReturnType
adj() const { return AdjReturnType(derived()); }


EIGEN_DEVICE_FUNC
inline NonConstAdjReturnType
adj() { return NonConstAdjReturnType(derived()); }