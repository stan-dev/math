#ifndef STAN_MATH_REV_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/trace_gen_quad_form.hpp>
#include <type_traits>

namespace stan {
namespace math {
namespace internal {

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
class trace_gen_quad_form_vari_alloc : public chainable_alloc {
 public:
  trace_gen_quad_form_vari_alloc(const Eigen::Matrix<Td, Rd, Cd>& D,
                                 const Eigen::Matrix<Ta, Ra, Ca>& A,
                                 const Eigen::Matrix<Tb, Rb, Cb>& B)
      : D_(D), A_(A), B_(B) {}

  double compute() {
    return trace_gen_quad_form(value_of(D_), value_of(A_), value_of(B_));
  }

  Eigen::Matrix<Td, Rd, Cd> D_;
  Eigen::Matrix<Ta, Ra, Ca> A_;
  Eigen::Matrix<Tb, Rb, Cb> B_;
};

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
class trace_gen_quad_form_vari : public vari {
 protected:
  static inline void computeAdjoints(double adj,
                                     const Eigen::Matrix<double, Rd, Cd>& D,
                                     const Eigen::Matrix<double, Ra, Ca>& A,
                                     const Eigen::Matrix<double, Rb, Cb>& B,
                                     Eigen::Matrix<var, Rd, Cd>* varD,
                                     Eigen::Matrix<var, Ra, Ca>* varA,
                                     Eigen::Matrix<var, Rb, Cb>* varB) {
    Eigen::Matrix<double, Ca, Cb> AtB;
    Eigen::Matrix<double, Ra, Cb> BD;
    if (varB || varA) {
      BD.noalias() = B * D;
    }
    if (varB || varD) {
      AtB.noalias() = A.transpose() * B;
    }

    if (varB) {
      (*varB).adj() += adj * (A * BD + AtB * D.transpose());
    }
    if (varA) {
      (*varA).adj() += adj * (B * BD.transpose());
    }
    if (varD) {
      (*varD).adj() += adj * (B.transpose() * AtB);
    }
  }

 public:
  explicit trace_gen_quad_form_vari(
      trace_gen_quad_form_vari_alloc<Td, Rd, Cd, Ta, Ra, Ca, Tb, Rb, Cb>* impl)
      : vari(impl->compute()), impl_(impl) {}

  virtual void chain() {
    computeAdjoints(adj_, value_of(impl_->D_), value_of(impl_->A_),
                    value_of(impl_->B_),
                    reinterpret_cast<Eigen::Matrix<var, Rd, Cd>*>(
                        std::is_same<Td, var>::value ? (&impl_->D_) : NULL),
                    reinterpret_cast<Eigen::Matrix<var, Ra, Ca>*>(
                        std::is_same<Ta, var>::value ? (&impl_->A_) : NULL),
                    reinterpret_cast<Eigen::Matrix<var, Rb, Cb>*>(
                        std::is_same<Tb, var>::value ? (&impl_->B_) : NULL));
  }

  trace_gen_quad_form_vari_alloc<Td, Rd, Cd, Ta, Ra, Ca, Tb, Rb, Cb>* impl_;
};
}  // namespace internal

template <typename Td, typename Ta, typename Tb,
          typename = require_any_var_t<value_type_t<Td>, value_type_t<Ta>,
                                       value_type_t<Tb>>,
          typename = require_all_eigen_t<Td, Ta, Tb>>
inline var trace_gen_quad_form(const Td& D, const Ta& A, const Tb& B) {
  using Td_scal = value_type_t<Td>;
  using Ta_scal = value_type_t<Ta>;
  using Tb_scal = value_type_t<Tb>;
  constexpr int Rd = Td::RowsAtCompileTime;
  constexpr int Cd = Td::ColsAtCompileTime;
  constexpr int Ra = Ta::RowsAtCompileTime;
  constexpr int Ca = Ta::ColsAtCompileTime;
  constexpr int Rb = Tb::RowsAtCompileTime;
  constexpr int Cb = Tb::ColsAtCompileTime;
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);

  auto* baseVari
      = new internal::trace_gen_quad_form_vari_alloc<Td_scal, Rd, Cd, Ta_scal,
                                                     Ra, Ca, Tb_scal, Rb, Cb>(
          D, A, B);

  return var(
      new internal::trace_gen_quad_form_vari<Td_scal, Rd, Cd, Ta_scal, Ra, Ca,
                                             Tb_scal, Rb, Cb>(baseVari));
}

template <typename Td, typename Ta, typename Tb,
          require_any_var_matrix_t<Td, Ta, Tb>* = nullptr,
          require_all_matrix_t<Td, Ta, Tb>* = nullptr>
inline var trace_gen_quad_form(const Td& D, const Ta& A, const Tb& B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);

  if(!is_constant<Ta>::value && !is_constant<Tb>::value && !is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    arena_t<promote_scalar_t<var, Ta>> arena_A = A;
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;

    var res = (arena_D.val_op() * arena_B.val_op().transpose() *
	       arena_A.val_op() * arena_B.val_op()).trace();

    reverse_pass_callback([arena_A, arena_B, arena_D, res]() mutable {
      double C_adj = res.adj();
      
      auto BDT = (arena_B.val_op() * arena_D.val_op().transpose()).eval();
      auto AB = (arena_A.val_op() * arena_B.val_op()).eval();

      arena_A.adj() += C_adj * BDT * arena_B.val_op().transpose();

      arena_B.adj()
	+= C_adj * (AB * arena_D.val_op() + arena_A.val_op().transpose() * BDT);
      
      arena_D.adj() += C_adj * (AB.transpose() * arena_B.val_op());
    });

    return res;
  } else if(!is_constant<Ta>::value && !is_constant<Tb>::value && is_constant<Td>::value) {
    arena_t<promote_scalar_t<var, Td>> arena_D = D;
    arena_t<promote_scalar_t<var, Ta>> arena_A = A;
    arena_t<promote_scalar_t<var, Tb>> arena_B = B;

    var res = (arena_D.val_op() * arena_B.val_op().transpose() *
	       arena_A.val_op() * arena_B.val_op()).trace();

    reverse_pass_callback([arena_A, arena_B, arena_D, res]() mutable {
      double C_adj = res.adj();
      
      auto BDT = (arena_B.val_op() * arena_D.val_op().transpose()).eval();

      arena_A.adj() += C_adj * BDT * arena_B.val_op().transpose();

      arena_B.adj()
	+= C_adj * (arena_A.val_op() * arena_B.val_op() * arena_D.val_op() +
		    arena_A.val_op().transpose() * BDT);
    });

    return res;
  } else if(!is_constant<Ta>::value && is_constant<Tb>::value && !is_constant<Td>::value) {
  } else if(!is_constant<Ta>::value && is_constant<Tb>::value && is_constant<Td>::value) {
  } else if(is_constant<Ta>::value && !is_constant<Tb>::value && !is_constant<Td>::value) {
  } else if(is_constant<Ta>::value && !is_constant<Tb>::value && is_constant<Td>::value) {
  } else if(is_constant<Ta>::value && is_constant<Tb>::value && !is_constant<Td>::value) {
  }
}

}  // namespace math
}  // namespace stan
#endif
