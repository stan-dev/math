/**
 * Structure to return value from an fvar. The first definition takes
 * a const fvar and returns a const double (for reading from a const matrix).
 * The second definition takes a non-const fvar and returns a non-const
 * double (for writing to a non-const matrix)
 */
struct val__Op {
  EIGEN_EMPTY_STRUCT_CTOR(val__Op)
  typedef typename Scalar::type result_type;
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const result_type& 
    operator()(const Scalar &v) const { return v.val_; }

  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE result_type& 
    operator()(Scalar &v) const { return v.val_; }
};

/**
 * Coefficient-wise function applying val__Op struct to a matrix of const fvar<T>
 * and returning a matrix of type T containing the values
 */
inline const CwiseUnaryView<val__Op, const Derived>
val_() const { return CwiseUnaryView<val__Op, const Derived>
    (derived(), val__Op());
}

/**
 * Coefficient-wise function applying val__Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the values that can
 * be modified
 */
inline CwiseUnaryView<val__Op, Derived>
val_() { return CwiseUnaryView<val__Op, Derived>
    (derived(), val__Op());
}

/**
 * Structure to return tangent from an fvar. The first definition takes
 * a const fvar and returns a const double (for reading from a const matrix).
 * The second definition takes a non-const fvar and returns a non-const
 * double (for writing to a non-const matrix)
 */
struct d__Op {
  EIGEN_EMPTY_STRUCT_CTOR(d__Op)
  typedef typename Scalar::type result_type;
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const result_type& 
    operator()(const Scalar &v) const { return v.d_; }

  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE result_type& 
    operator()(Scalar &v) const { return v.d_; }
};

/**
 * Coefficient-wise function applying d__Op struct to a matrix of const fvar<T>
 * and returning a const matrix of type T containing the tangents
 */
inline const CwiseUnaryView<d__Op, const Derived>
d_() const { return CwiseUnaryView<d__Op, const Derived>
    (derived(), d__Op());
}

/**
 * Coefficient-wise function applying d__Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the tangents that can
 * be modified
 */
inline CwiseUnaryView<d__Op, Derived>
d_() { return CwiseUnaryView<d__Op, Derived>
    (derived(), d__Op());
}

/**
 * Structure to return value from a var. The definition takes
 * a const var and returns a const double (for reading from a const matrix).
 */
struct val_Op {
  EIGEN_EMPTY_STRUCT_CTOR(val_Op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const double& 
    operator()(const Scalar &v) const { return v.vi_->val_; }
};

/**
 * Coefficient-wise function applying val_Op struct to a matrix of const var
 * and returning a const matrix of type T containing the values
 */
inline const CwiseUnaryView<val_Op, const Derived>
val() const { return CwiseUnaryView<val_Op, const Derived>
    (derived(), val_Op());
}

/**
 * Structure to return adjoint from a var. The first definition takes
 * a const var and returns a const double (for reading from a const matrix).
 * The second definition takes a non-const var and returns a non-const
 * double (for writing to a non-const matrix)
 */
struct adj_Op {
  EIGEN_EMPTY_STRUCT_CTOR(adj_Op)
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const double& 
    operator()(const Scalar &v) const { return v.vi_->adj_; }

  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE double& 
    operator()(Scalar &v) const { return v.vi_->adj_; }
};

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of const var
 * and returning a const matrix of type T containing the values
 */
inline const CwiseUnaryView<adj_Op, const Derived>
adj() const { return CwiseUnaryView<adj_Op, const Derived>
    (derived(), adj_Op());
}

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of var
 * and returning a view to a matrix of type T of the adjoints that can
 * be modified
 */
inline CwiseUnaryView<adj_Op, Derived>
adj() { return CwiseUnaryView<adj_Op, Derived>
    (derived(), adj_Op());
}
