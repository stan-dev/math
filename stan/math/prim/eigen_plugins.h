// Workaround for SFINAE regression with GCC 8.3
template <typename> struct has_member { typedef int type; };

/**
 * Reimplements is_fvar without requiring external math headers
 *
 * decltype((void)(T::d_)) is a pre C++17 replacement for
 * std::void_t<decltype(T::d_)>
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template <typename T,
          typename has_member<decltype(std::decay_t<T>::d_)>::type = 0>
struct is_fvar : std::true_type
{ };

/**
 * Reimplements is_var without requiring external math headers
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template <typename T,
          typename has_member<decltype(std::decay_t<T>::vi_)>::type = 0>
struct is_var : std::true_type
{ };

/**
 * Reimplements is_vari without requiring external math headers
 *
 * TODO(Andrew): Replace with std::void_t after move to C++17
 */
template <typename T,
          typename has_member<
                      decltype(std::remove_pointer_t<
                                  std::decay_t<T>>::adj_)>::type = 0>
struct is_vari : std::true_type
{ };


/**
 * TODO(Andrew): Replace std::is_const<>::value with std::is_const_v<>
 *                after move to C++17
 */
template <typename T_in, typename T_return>
using match_const_t =
        std::conditional_t<std::is_const<std::remove_reference_t<T_in>>::value,
                           const T_return&, T_return&>;

template <typename T>
using reverse_return_t = match_const_t<T, double>;

template <typename T>
using vari_return_t = match_const_t<T, decltype(std::decay_t<T>::vi_)>;

template <typename T>
using forward_return_t = match_const_t<T, decltype(std::decay_t<T>::val_)>;

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
  template <typename T = Scalar,
           std::enable_if_t<is_vari<T>::value>* = nullptr>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    reverse_return_t<T> operator()(T&& v) const { return v->val_; }

  //Returns value from a var
  template <typename T = Scalar,
           std::enable_if_t<is_var<T>::value>* = nullptr>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    reverse_return_t<T> operator()(T&& v) const {
      return v.vi_->val_;
    }

  //Returns value from an fvar
  template <typename T = Scalar,
           std::enable_if_t<is_fvar<T>::value>* = nullptr>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    forward_return_t<T> operator()(T&& v) const { return v.val_; }

  //Returns double unchanged from input (by reference)
  template <typename T = Scalar,
           std::enable_if_t<std::is_arithmetic<std::decay_t<T>>::value>* = nullptr>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  match_const_t<T, double> operator()(T&& v) const { return v; }
};

/**
 * Coefficient-wise function applying val_Op struct to a matrix of const var,
 * vari*, fvar, or doubles and returning a view to the const matrix
 * containing the values
 */
inline const CwiseUnaryOp<val_Op, Derived>
val() const { return CwiseUnaryOp<val_Op, Derived>(derived()); }

/**
 * Coefficient-wise function applying val_Op struct to a matrix of var
 * or vari* and returning a view to the const matrix containing
 * the values
 */
inline const CwiseUnaryOp<val_Op, Derived>
val_op() const { return CwiseUnaryOp<val_Op, Derived>(derived()); }

/**
 * Coefficient-wise function applying val_Op struct to a matrix of var
 * or vari* and returning a view to the values
 */
inline CwiseUnaryView<val_Op, Derived>
val() { return CwiseUnaryView<val_Op, Derived>(derived()); }

/**
 * Structure to return tangent from an fvar.
 */
struct d_Op {
  EIGEN_EMPTY_STRUCT_CTOR(d_Op);

  //Returns tangent from an fvar
  template <typename T = Scalar,
            std::enable_if_t<is_fvar<T>::value>* = nullptr>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    forward_return_t<T> operator()(T&& v) const { return v.d_; }
};

/**
 * Coefficient-wise function applying d_Op struct to a matrix of const fvar<T>
 * and returning a const matrix of type T containing the tangents
 */
inline const CwiseUnaryOp<d_Op, Derived>
d() const { return CwiseUnaryOp<d_Op, Derived>(derived()); }

/**
 * Coefficient-wise function applying d_Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the tangents.
 */
inline const CwiseUnaryOp<d_Op, Derived>
d_op() const { return CwiseUnaryOp<d_Op, Derived>(derived()); }

/**
 * Coefficient-wise function applying d_Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the tangents that can
 * be modified
 */
inline CwiseUnaryView<d_Op, Derived>
d() { return CwiseUnaryView<d_Op, Derived>(derived()); }

/**
 * Structure to return adjoints from var and vari*.
 */
struct adj_Op {
  EIGEN_EMPTY_STRUCT_CTOR(adj_Op);

  //Returns adjoint from a vari*
  template <typename T = Scalar,
           std::enable_if_t<is_vari<T>::value>* = nullptr>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    reverse_return_t<T> operator()(T&& v) const { return v->adj_; }

  //Returns adjoint from a var
  template <typename T = Scalar,
           std::enable_if_t<is_var<T>::value>* = nullptr>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    reverse_return_t<T> operator()(T&& v) const { return v.vi_->adj_; }
};

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of const var
 * and returning a const matrix of type T containing the values
 */
inline const CwiseUnaryOp<adj_Op, Derived>
adj() const { return CwiseUnaryOp<adj_Op, Derived>(derived()); }

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of var
 * and returning a view to a matrix of doubles of the adjoints.
 * This is meant to be used on the rhs of expressions.
 */
inline const CwiseUnaryOp<adj_Op, Derived>
adj_op() const { return CwiseUnaryOp<adj_Op, Derived>(derived()); }

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of var
 * and returning a view to a matrix of doubles of the adjoints that can
 * be modified
 */
inline CwiseUnaryView<adj_Op, Derived>
adj() { return CwiseUnaryView<adj_Op, Derived>(derived()); }

/**
 * Structure to return vari* from a var.
 */
struct vi_Op {
  EIGEN_EMPTY_STRUCT_CTOR(vi_Op);

  //Returns vari* from a var
  template <typename T = Scalar>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    vari_return_t<T> operator()(T&& v) const { return v.vi_; }
};

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of const var
 * and returning a const matrix of vari*
 */
inline const CwiseUnaryOp<vi_Op, Derived>
vi() const { return CwiseUnaryOp<vi_Op, Derived>(derived()); }

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of var
 * and returning a view to a matrix of vari* that can be modified
 */
inline CwiseUnaryView<vi_Op, Derived>
vi() { return CwiseUnaryView<vi_Op, Derived>(derived()); }

#define EIGEN_STAN_MATRIXBASE_PLUGIN
#define EIGEN_STAN_ARRAYBASE_PLUGIN
