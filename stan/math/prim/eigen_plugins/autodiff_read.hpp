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

template <typename EigVari, typename EigDbl>
inline void read_vi_val_adj(EigVari& VariMat, EigDbl& ValMat,
                            EigDbl& AdjMat) const {
    AdjMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      vi_val_adj_functor<EigVari, EigDbl>(this->derived(), VariMat, ValMat)
    );
}

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

template <typename EigDbl>
inline void read_val_adj(EigDbl& ValMat, EigDbl& AdjMat) const {
    AdjMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      val_adj_functor<EigDbl>(derived(), ValMat)
    );
}

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
    auto operator() (Index row, Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    return var_mat.coeffRef(row, col).vi_->val_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    auto operator() (Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    return var_mat.coeffRef(index).vi_->val_;
  }
};

template <typename EigVari, typename EigDbl>
inline void read_vi_val(EigVari& VariMat, EigDbl& ValMat) const {
    ValMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      vi_val_functor<EigVari>(this->derived(), VariMat)
    );
}

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
    auto operator() (Index row, Index col) const {
    vi_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).vi_;
    return var_mat.coeffRef(row, col).vi_->adj_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    auto operator() (Index index) const {
    vi_mat.coeffRef(index) = var_mat.coeffRef(index).vi_;
    return var_mat.coeffRef(index).vi_->adj_;
  }
};

template <typename EigVari, typename EigDbl>
inline void read_vi_adj(EigVari& VariMat, EigDbl& AdjMat) const {
    AdjMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      vi_val_functor<EigVari>(this->derived(), VariMat)
    );
}

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
    auto operator() (Index row, Index col) const {
    val_mat.coeffRef(row, col) = var_mat.coeffRef(row, col).val_;
    return var_mat.coeffRef(row, col).d_;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
    auto operator() (Index index) const {
    val_mat.coeffRef(index) = var_mat.coeffRef(index).val_;
    return var_mat.coeffRef(index).d_;
  }
};

template <typename EigDbl>
inline void read_fvar(EigDbl& ValMat, EigDbl& DMat) const {
    DMat = EigDbl::NullaryExpr(
      this->rows(),
      this->cols(),
      read_fvar_functor<EigDbl>(this->derived(), ValMat)
    );
}