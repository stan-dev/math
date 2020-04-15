/**
 * Reimplements is_fvar without requiring external math headers
 *
 * decltype((void)(T::d_)) is a pre C++17 replacement for
 * std::void_t<decltype(T::d_)>
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template<class, class = void>
struct is_fvar : std::false_type
{ };
template<class T>
struct is_fvar<T, decltype((void)(T::d_))> : std::true_type
{ };

/**
 * Reimplements is_var without requiring external math headers
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template<class, class = void>
struct is_var : std::false_type
{ };
template<class T>
struct is_var<T, decltype((void)(T::vi_))> : std::true_type
{ };

//TODO(Andrew): Replace std::is_const<>::value with std::is_const_v<> after move to C++17
template<typename T>
using double_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                         const double,
                                         double>;
template<typename T>
using reverse_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                         const double&,
                                         double&>;

template<typename T>
using vari_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                          const decltype(T::vi_)&,
                                          decltype(T::vi_)&>;

template<typename T>
using forward_return_t = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                                         const decltype(T::val_)&,
                                         decltype(T::val_)&>;

template <typename T>
using require_double_t = std::enable_if_t<std::is_floating_point<T>::value,double_return_t<T>>;

template <typename T>
using require_var_t = std::enable_if_t<is_var<T>::value,reverse_return_t<T>>;

template <typename T>
using require_fvar_t = std::enable_if_t<is_fvar<T>::value,forward_return_t<T>>;

template <typename T>
using require_vari_t = std::enable_if_t<std::is_pointer<T>::value,reverse_return_t<T>>;

/**
 * Structure to return a view to the values in a var, vari*, and fvar<T>.
 * To identify the correct member to call for a given input, templates
 * check a combination of whether the input is a pointer (i.e. vari*)
 * and/or whether the input has member ".d_" (i.e. fvar).
 *
 * There are two methods for returning doubles unchanged. One which takes a reference
 * to a double and returns the same reference, used when 'chaining' methods
 * (i.e. A.adj().val()). The other for passing and returning by value, used directly
 * with matrices of doubles (i.e. A.val(), where A is of type MatrixXd).
 *
 * For definitions of EIGEN_EMPTY_STRUCT_CTOR, EIGEN_DEVICE_FUNC, and
 * EIGEN_STRONG_INLINE; see: https://eigen.tuxfamily.org/dox/XprHelper_8h_source.html
 */
struct val_Op{
  EIGEN_EMPTY_STRUCT_CTOR(val_Op);

  //Returns value from a vari*
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    require_vari_t<T> operator()(T& v) const { return v->val_; }

  //Returns value from a var
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    require_var_t<T> operator()(T& v) const { return v.vi_->val_; }

  //Returns value from an fvar
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  require_fvar_t<T> operator()(T& v) const { return v.val_; }

  //Returns double unchanged from input (by value)
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    require_double_t<T> operator()(T v) const { return v; }

  //Returns double unchanged from input (by reference)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  const double& operator()(const double& v) const { return v; }

  //Returns double unchanged from input (by reference)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  double& operator()(double& v) const { return v; }
};

/**
 * Coefficient-wise function applying val_Op struct to a matrix of const var
 * or vari* and returning a view to the const matrix of doubles containing
 * the values
 */
inline const CwiseUnaryOp<val_Op, const Derived>
val() const { return CwiseUnaryOp<val_Op, const Derived>(derived());
}

/**
 * Coefficient-wise function applying val_Op struct to a matrix of var
 * or vari* and returning a view to the values
 */
inline CwiseUnaryView<val_Op, Derived>
val() { return CwiseUnaryView<val_Op, Derived>(derived());
}

/**
 * Structure to return tangent from an fvar.
 */
struct d_Op {
  EIGEN_EMPTY_STRUCT_CTOR(d_Op);

  //Returns tangent from an fvar
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    forward_return_t<T> operator()(T& v) const { return v.d_; }
};

/**
 * Coefficient-wise function applying d_Op struct to a matrix of const fvar<T>
 * and returning a const matrix of type T containing the tangents
 */
inline const CwiseUnaryOp<d_Op, const Derived>
d() const { return CwiseUnaryOp<d_Op, const Derived>(derived());
}

/**
 * Coefficient-wise function applying d_Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the tangents that can
 * be modified
 */
inline CwiseUnaryView<d_Op, Derived>
d() { return CwiseUnaryView<d_Op, Derived>(derived());
}

/**
 * Structure to return adjoints from var and vari*. Deduces whether the variables
 * are pointers (i.e. vari*) to determine whether to return the adjoint or
 * first point to the underlying vari* (in the case of var).
 */
struct adj_Op {
  EIGEN_EMPTY_STRUCT_CTOR(adj_Op);

  //Returns adjoint from a vari*
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    require_vari_t<T> operator()(T& v) const { return v->adj_; }

  //Returns adjoint from a var
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    require_var_t<T> operator()(T& v) const { return v.vi_->adj_; }
};

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of const var
 * and returning a const matrix of type T containing the values
 */
inline const CwiseUnaryOp<adj_Op, const Derived>
adj() const { return CwiseUnaryOp<adj_Op, const Derived>(derived());
}

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of var
 * and returning a view to a matrix of doubles of the adjoints that can
 * be modified
 */
inline CwiseUnaryView<adj_Op, Derived>
adj() { return CwiseUnaryView<adj_Op, Derived>(derived());
}
/**
 * Structure to return vari* from a var.
 */
struct vi_Op {
  EIGEN_EMPTY_STRUCT_CTOR(vi_Op);

  //Returns vari* from a var
  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    vari_return_t<T> operator()(T& v) const { return v.vi_; }
};

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of const var
 * and returning a const matrix of vari*
 */
inline const CwiseUnaryOp<vi_Op, const Derived>
vi() const { return CwiseUnaryOp<vi_Op, const Derived>(derived());
}

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of var
 * and returning a view to a matrix of vari* that can be modified
 */
inline CwiseUnaryView<vi_Op, Derived>
vi() { return CwiseUnaryView<vi_Op, Derived>(derived());
}

/**
 * Functor for extracting the vari*, values, and adjoints from a matrix of var.
 * This functor is called using Eigen's NullaryExpr framework, which takes care
 * of the indexing. This removes the need to programmatically account for
 * whether the input is row- or column-major.
 */
template<typename EigVari, typename EigDbl>
class vi_val_adj_functor
{
  const Derived& var_mat;
  EigVari& vi_mat;
  EigDbl& val_mat;

public:
  vi_val_adj_functor(const Derived& arg1, EigVari& arg2, EigDbl& arg3)
    : var_mat(arg1), vi_mat(arg2), val_mat(arg3)
  {}

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    auto operator() (Index row, Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_->val_;
    return var_mat.coeffRef(row, col).vi_->adj_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    auto operator() (Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    val_mat.coeffRef(index) = var_mat.coeffRef(index).vi_->val_;
    return var_mat.coeffRef(index).vi_->adj_;
  }
};

/**
 * Member function applying the vi_val_adj_functor to extract the vari*, values,
 * and adjoints of a given var matrix into separate matrices.
 */
template <typename EigVari, typename EigDbl>
inline void read_vi_val_adj(EigVari& VariMat, EigDbl& ValMat,
                            EigDbl& AdjMat) const {
    AdjMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      vi_val_adj_functor<EigVari, EigDbl>(this->derived(), VariMat, ValMat)
    );
}

/**
 * Functor for extracting the values and adjoints from a matrix of var or vari.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template<typename EigDbl>
class val_adj_functor
{
  const Derived& var_mat;
  EigDbl& val_mat;

public:
  val_adj_functor(const Derived& arg1, EigDbl& arg2)
    : var_mat(arg1), val_mat(arg2)
  {}

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  require_var_t<T> operator() (Index row, Index col) const {
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_->val_;
    return var_mat.coeffRef(row, col).vi_->adj_;
  }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  require_var_t<T> operator() (Index index) const {
    val_mat.coeffRef(index) = var_mat.coeffRef(index).vi_->val_;
    return var_mat.coeffRef(index).vi_->adj_;
  }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  require_vari_t<T> operator() (Index row, Index col) const {
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).val_;
    return var_mat.coeffRef(row, col).adj_;
  }

  template<typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  require_vari_t<T> operator() (Index index) const {
    val_mat.coeffRef(index) = var_mat.coeffRef(index).val_;
    return var_mat.coeffRef(index).adj_;
  }
};

/**
 * Member function applying the val_adj_functor to extract the values
 * and adjoints of a given var or vari matrix into separate matrices.
 */
template <typename EigDbl>
inline void read_val_adj(EigDbl& ValMat, EigDbl& AdjMat) const {
    AdjMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      val_adj_functor<EigDbl>(derived(), ValMat)
    );
}

/**
 * Functor for extracting the varis and values from a matrix of var.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template<typename EigVari>
class vi_val_functor
{
  const Derived& var_mat;
  EigVari& vi_mat;

public:
  vi_val_functor(const Derived& arg1, EigVari& arg2)
    : var_mat(arg1), vi_mat(arg2)
  {}

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    decltype(auto) operator() (Index row, Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    return var_mat.coeffRef(row, col).vi_->val_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    decltype(auto) operator() (Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    return var_mat.coeffRef(index).vi_->val_;
  }
};

/**
 * Member function applying the vi_val_functor to extract the varis
 * and values of a given var matrix into separate matrices.
 */
template <typename EigVari, typename EigDbl>
inline void read_vi_val(EigVari& VariMat, EigDbl& ValMat) const {
    ValMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      vi_val_functor<EigVari>(this->derived(), VariMat)
    );
}

/**
 * Functor for extracting the varis and adjoints from a matrix of var.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template<typename EigVari>
class vi_adj_functor
{
  const Derived& var_mat;
  EigVari& vi_mat;

public:
  vi_adj_functor(const Derived& arg1, EigVari& arg2)
    : var_mat(arg1), vi_mat(arg2)
  {}

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    decltype(auto) operator() (Index row, Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    return var_mat.coeffRef(row, col).vi_->adj_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    decltype(auto) operator() (Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    return var_mat.coeffRef(index).vi_->adj_;
  }
};

/**
 * Member function applying the vi_adj_functor to extract the varis
 * and adjoints of a given var matrix into separate matrices.
 */
template <typename EigVari, typename EigDbl>
inline void read_vi_adj(EigVari& VariMat, EigDbl& AdjMat) const {
    AdjMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      vi_val_functor<EigVari>(this->derived(), VariMat)
    );
}

/**
 * Functor for extracting the values and tangents from a matrix of fvar.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template<typename EigDbl>
class read_fvar_functor
{
  const Derived& var_mat;
  EigDbl& val_mat;

public:
  read_fvar_functor(const Derived& arg1, EigDbl& arg2)
    : var_mat(arg1), val_mat(arg2)
  {}

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    decltype(auto) operator() (Index row, Index col) const {
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).val_;
    return var_mat.coeffRef(row, col).d_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    decltype(auto) operator() (Index index) const {
    val_mat.coeffRef(index) = var_mat.coeffRef(index).val_;
    return var_mat.coeffRef(index).d_;
  }
};

/**
 * Member function applying the read_fvar_functor to extract the values
 * and tangents of a given fvar matrix into separate matrices.
 */
template <typename EigDbl>
inline void read_fvar(EigDbl& ValMat, EigDbl& DMat) const {
    DMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      read_fvar_functor<EigDbl>(this->derived(), ValMat)
    );
}

#define EIGEN_STAN_MATRIXBASE_PLUGIN
#define EIGEN_STAN_ARRAYBASE_PLUGIN
