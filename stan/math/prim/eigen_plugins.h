

/**
 * Coefficient-wise function applying val_Op struct to a matrix of const var,
 * vari*, fvar, or doubles and returning a view to the const matrix
 * containing the values
 */
EIGEN_DEVICE_FUNC inline const auto val() const {
   return CwiseUnaryOp<const internal::val_Op, const Derived>(derived());
}

/**
 * Coefficient-wise function applying val_Op struct to a matrix of var
 * or vari* and returning a view to the values
 */
EIGEN_DEVICE_FUNC inline auto val() {
   return CwiseUnaryView<internal::val_ref_Op, Derived>(derived());
}

/**
 * Coefficient-wise function applying d_Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the tangents that can
 * be modified
 */
EIGEN_DEVICE_FUNC inline auto d() const {
  return CwiseUnaryView<internal::d_ref_Op, Derived>(const_cast_derived());
}

/**
 * Coefficient-wise function applying d_Op struct to a matrix of fvar<T>
 * and returning a view to a matrix of type T of the tangents that can
 * be modified
 */
EIGEN_DEVICE_FUNC inline auto d() {
  return CwiseUnaryView<internal::d_ref_Op, Derived>(derived());
}

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of const var
 * and returning a const matrix of type T containing the values
 */
EIGEN_DEVICE_FUNC inline const auto adj() const {
  return CwiseUnaryOp<const internal::adj_Op, const Derived>(derived());
}

/**
 * Coefficient-wise function applying adj_Op struct to a matrix of var
 * and returning a view to a matrix of doubles of the adjoints that can
 * be modified
 */
EIGEN_DEVICE_FUNC inline auto adj() {
   return CwiseUnaryView<internal::adj_ref_Op, Derived>(derived());
}

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of const var
 * and returning a const matrix of vari*
 */
EIGEN_DEVICE_FUNC inline const auto vi() const {
  return CwiseUnaryOp<const internal::vi_Op, const Derived>(derived());
}

/**
 * Coefficient-wise function applying vi_Op struct to a matrix of var
 * and returning a view to a matrix of vari* that can be modified
 */
EIGEN_DEVICE_FUNC inline auto vi() {
  return CwiseUnaryView<internal::vi_ref_Op, Derived>(derived());
}

#define EIGEN_STAN_MATRIXBASE_PLUGIN
#define EIGEN_STAN_ARRAYBASE_PLUGIN
