#ifndef STAN_MATH_REV_CORE_READ_VAR_HPP
#define STAN_MATH_REV_CORE_READ_VAR_HPP

#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>

namespace stan {
namespace math {

/**
 * Functor for extracting the vari*, values, and adjoints from a matrix of var.
 * This functor is called using Eigen's NullaryExpr framework, which takes care
 * of the indexing. This removes the need to programmatically account for
 * whether the input is row- or column-major.
 */
template <typename EigRev, typename EigVari, typename EigDbl>
class vi_val_adj_functor {
  const EigRev& var_mat;
  EigVari& vi_mat;
  EigDbl& val_mat;

 public:
  vi_val_adj_functor(const EigRev& arg1, EigVari& arg2, EigDbl& arg3)
      : var_mat(arg1), vi_mat(arg2), val_mat(arg3) {}

  inline decltype(auto) operator()(Eigen::Index row, Eigen::Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_->val_;
    return var_mat.coeffRef(row, col).vi_->adj_;
  }

  inline decltype(auto) operator()(Eigen::Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    val_mat.coeffRef(index) = var_mat.coeffRef(index).vi_->val_;
    return var_mat.coeffRef(index).vi_->adj_;
  }
};

/**
 * Functor for extracting the values and adjoints from a matrix of var or vari.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template <typename EigRev, typename EigDbl>
class val_adj_functor {
  const EigRev& var_mat;
  EigDbl& val_mat;

 public:
  val_adj_functor(const EigRev& arg1, EigDbl& arg2)
      : var_mat(arg1), val_mat(arg2) {}

  template <typename T = EigRev, require_st_same<T, var>* = nullptr>
  inline decltype(auto) operator()(Eigen::Index row, Eigen::Index col) const {
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_->val_;
    return var_mat.coeffRef(row, col).vi_->adj_;
  }

  template <typename T = EigRev, require_st_same<T, var>* = nullptr>
  inline decltype(auto) operator()(Eigen::Index index) const {
    val_mat.coeffRef(index) = var_mat.coeffRef(index).vi_->val_;
    return var_mat.coeffRef(index).vi_->adj_;
  }

  template <typename T = EigRev, require_st_same<T, vari*>* = nullptr>
  inline decltype(auto) operator()(Eigen::Index row, Eigen::Index col) const {
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col)->val_;
    return var_mat.coeffRef(row, col)->adj_;
  }

  template <typename T = EigRev, require_st_same<T, vari*>* = nullptr>
  inline decltype(auto) operator()(Eigen::Index index) const {
    val_mat.coeffRef(index) = var_mat.coeffRef(index)->val_;
    return var_mat.coeffRef(index)->adj_;
  }
};

/**
 * Functor for extracting the varis and values from a matrix of var.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template <typename EigVar, typename EigVari>
class vi_val_functor {
  const EigVar& var_mat;
  EigVari& vi_mat;

 public:
  vi_val_functor(const EigVar& arg1, EigVari& arg2)
      : var_mat(arg1), vi_mat(arg2) {}

  inline decltype(auto) operator()(Eigen::Index row, Eigen::Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    return var_mat.coeffRef(row, col).vi_->val_;
  }

  inline decltype(auto) operator()(Eigen::Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    return var_mat.coeffRef(index).vi_->val_;
  }
};

/**
 * Functor for extracting the varis and adjoints from a matrix of var.
 * This functor is called using Eigen's NullaryExpr framework.
 */
template <typename EigVar, typename EigVari>
class vi_adj_functor {
  const EigVar& var_mat;
  EigVari& vi_mat;

 public:
  vi_adj_functor(const EigVar& arg1, EigVari& arg2)
      : var_mat(arg1), vi_mat(arg2) {}

  inline decltype(auto) operator()(Eigen::Index row, Eigen::Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    return var_mat.coeffRef(row, col).vi_->adj_;
  }

  inline decltype(auto) operator()(Eigen::Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    return var_mat.coeffRef(index).vi_->adj_;
  }
};

/**
 * Function applying the vi_val_adj_functor to extract the vari*, values,
 * and adjoints of a given var matrix into separate matrices.
 *
 * @tparam EigVar type of the Eigen container of var.
 * @tparam EigVari type of the Eigen container of vari to be copied to.
 * @tparam EigDbl type of the Eigen container of doubles to be copied to.
 * @param[in] VarMat Input Eigen container of var.
 * @param[in] VariMat Output Eigen container of vari.
 * @param[in] ValMat Output Eigen container of values.
 * @param[in] AdjMat Output Eigen container of tangents.
 */
template <typename EigVar, typename EigVari, typename EigDbl>
inline void read_vi_val_adj(const EigVar& VarMat, EigVari& VariMat,
                            EigDbl& ValMat, EigDbl& AdjMat) {
  AdjMat
      = EigDbl::NullaryExpr(VarMat.rows(), VarMat.cols(),
                            vi_val_adj_functor<const EigVar, EigVari, EigDbl>(
                                VarMat, VariMat, ValMat));
}

/**
 * Function applying the val_adj_functor to extract the values
 * and adjoints of a given var or vari matrix into separate matrices.
 *
 * @tparam EigRev type of the Eigen container of var or vari.
 * @tparam EigDbl type of the Eigen container of doubles to be copied to.
 * @param[in] VarMat Input Eigen container of var.
 * @param[in] ValMat Output Eigen container of values.
 * @param[in] AdjMat Output Eigen container of adjoints.
 */
template <typename EigRev, typename EigDbl>
inline void read_val_adj(const EigRev& VarMat, EigDbl& ValMat, EigDbl& AdjMat) {
  AdjMat = EigDbl::NullaryExpr(
      VarMat.rows(), VarMat.cols(),
      val_adj_functor<const EigRev, EigDbl>(VarMat, ValMat));
}

/**
 * Function applying the vi_val_functor to extract the varis and
 * and values of a given var matrix into separate matrices.
 *
 * @tparam EigVar type of the Eigen container of var.
 * @tparam EigDbl type of the Eigen container of doubles to be copied to.
 * @param[in] VarMat Input Eigen container of var.
 * @param[in] VariMat Output Eigen container of vari.
 * @param[in] ValMat Output Eigen container of values.
 */
template <typename EigVar, typename EigVari, typename EigDbl>
inline void read_vi_val(const EigVar& VarMat, EigVari& VariMat,
                        EigDbl& ValMat) {
  ValMat = EigDbl::NullaryExpr(
      VarMat.rows(), VarMat.cols(),
      vi_val_functor<const EigVar, EigVari>(VarMat, VariMat));
}

/**
 * Function applying the vi_adj_functor to extract the varis and
 * and adjoints of a given var matrix into separate matrices.
 *
 * @tparam EigVar type of the Eigen container of var.
 * @tparam EigDbl type of the Eigen container of doubles to be copied to.
 * @param[in] VarMat Input Eigen container of var.
 * @param[in] VariMat Output Eigen container of vari.
 * @param[in] AdjMat Output Eigen container of adjoints.
 */
template <typename EigVar, typename EigVari, typename EigDbl>
inline void read_vi_adj(const EigVar& VarMat, EigVari& VariMat,
                        EigDbl& AdjMat) {
  AdjMat = EigDbl::NullaryExpr(
      VarMat.rows(), VarMat.cols(),
      vi_adj_functor<const EigVar, EigVari>(VarMat, VariMat));
}

}  // namespace math
}  // namespace stan
#endif
