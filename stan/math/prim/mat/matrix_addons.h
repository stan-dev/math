/**
 * Structure to determine whether input struct has member named "d_".
 * This is used to differentiate between var and fvar<T> without
 * depending on external code.
 *
 * Re-implements std::void_t for pre-c++17 to detect ill-formed types
 * in SFINAE
 */
template<class...>
using void_t = void;

template<class, class = void>
struct Hasd_ : std::false_type
{ };
template<class T>
struct Hasd_<T, void_t<decltype(T::d_)>> : std::true_type
{ };

/**
 * Structure to return a view to the values in a var, vari*, and fvar<T>.
 * To identify the correct member to call for a given input, templates
 * check a combination of whether the input is a pointer (i.e. vari*)
 * and/or whether the input has member ".d_" (i.e. fvar).
 *
 * The operators are overloaded for both const and non-const inputs
 */
struct val_Op{
  EIGEN_EMPTY_STRUCT_CTOR(val_Op);
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename std::enable_if_t<std::is_pointer<T>::value, const double&>
      operator()(const Scalar &v) const { return v->val_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename std::enable_if_t<(!std::is_pointer<T>::value && !Hasd_<T>::value),
                               const double&>
      operator()(const Scalar &v) const { return v.vi_->val_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename std::enable_if_t<Hasd_<T>::value, const decltype(T::val_)&>
      operator()(const Scalar &v) const { return v.val_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename std::enable_if_t<std::is_pointer<T>::value, double&>
      operator()(Scalar &v) const { return v->val_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename std::enable_if_t<(!std::is_pointer<T>::value && !Hasd_<T>::value),
                                                double&>
      operator()(Scalar &v) { return v.vi_->val_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    typename std::enable_if_t<Hasd_<T>::value, decltype(T::val_)&>
      operator()(Scalar &v) { return v.val_; }
};

/**
 * Coefficient-wise function applying val_Op struct to a matrix of const var
 * or vari* and returning a view to the const matrix of doubles containing
 * the values
 */
inline const CwiseUnaryOp<val_Op, const Derived>
val() const { return CwiseUnaryOp<val_Op, const Derived>
    (derived(), val_Op());
}

/**
 * Coefficient-wise function applying val_Op struct to a matrix of const var
 * or vari* and returning a view to the matrix of doubles containing the
 * values
 */
inline CwiseUnaryView<val_Op, Derived>
val() { return CwiseUnaryView<val_Op, Derived>
    (derived(), val_Op());
}

/**
 * Structure to return tangent from an fvar. The first definition takes
 * a const fvar and returns a const double (for reading from a const matrix).
 * The second definition takes a non-const fvar and returns a non-const
 * double (for writing to a non-const matrix)
 */
struct d_Op {
  EIGEN_EMPTY_STRUCT_CTOR(d_Op);
  typedef decltype(Scalar::d_) result_type;
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
inline const CwiseUnaryOp<d_Op, const Derived>
d() const { return CwiseUnaryOp<d_Op, const Derived>
    (derived(), d_Op());
}

/**
 * Coefficient-wise function applying d__Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the tangents that can
 * be modified
 */
inline CwiseUnaryView<d_Op, Derived>
d() { return CwiseUnaryView<d_Op, Derived>
    (derived(), d_Op());
}

/**
 * Structure to return adjoints from var and vari*. Tests whether the variables
 * are pointers (i.e. vari*) to determine whether to return the adjoint or
 * first point to the underlying vari* (in the case of var). The first
 * definition takes a const var and returns a const double (for reading from
 * a const matrix). The second definition takes a non-const var and returns
 * a non-const double (for writing to a non-const matrix)
 */
struct adj_Op {
  EIGEN_EMPTY_STRUCT_CTOR(adj_Op);
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE typename
    std::enable_if_t<std::is_pointer<T>::value,const double&> const 
    operator()(const Scalar &v) const { return v->adj_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE typename
    std::enable_if_t<!std::is_pointer<T>::value,const double&> const 
    operator()(const Scalar &v) const { return v.vi_->adj_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE typename
    std::enable_if_t<std::is_pointer<T>::value,double&>
    operator()(Scalar &v) const { return v->adj_; }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE typename
    std::enable_if_t<!std::is_pointer<T>::value,double&>
    operator()(Scalar &v) const { return v.vi_->adj_; }
};

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of const var
 * and returning a const matrix of type T containing the values
 */
inline const CwiseUnaryOp<adj_Op, const Derived>
adj() const { return CwiseUnaryOp<adj_Op, const Derived>
    (derived(), adj_Op());
}

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of var
 * and returning a view to a matrix of doubles of the adjoints that can
 * be modified
 */
inline CwiseUnaryView<adj_Op, Derived>
adj() { return CwiseUnaryView<adj_Op, Derived>
    (derived(), adj_Op());
}
/**
 * Structure to return vari* from a var. The first definition takes
 * a const var and returns a const vari* (for reading from a const matrix).
 * The second definition takes a non-const var and returns a non-const
 * vari* (for writing to a non-const matrix)
 */
struct vi_Op {
  typedef decltype(Scalar::vi_) result_type;
  EIGEN_EMPTY_STRUCT_CTOR(vi_Op);
  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE const result_type& 
    operator()(const Scalar &v) const { return v.vi_; }

  EIGEN_DEVICE_FUNC
  EIGEN_STRONG_INLINE result_type& 
    operator()(Scalar &v) const { return v.vi_; }
};

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of const var
 * and returning a const matrix of vari*
 */
inline const CwiseUnaryOp<vi_Op, const Derived>
vi() const { return CwiseUnaryOp<vi_Op, const Derived>
    (derived(), vi_Op());
}

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of var
 * and returning a view to a matrix of vari* that can be modified
 */
inline CwiseUnaryView<vi_Op, Derived>
vi() { return CwiseUnaryView<vi_Op, Derived>
    (derived(), vi_Op());
}
