#ifndef STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <iostream>
#include <iterator>
#include <vector>

namespace stan {
namespace math {

template <typename T>
struct value_type {
  using type = double;
};

template <>
struct value_type<int> {
  using type = int;
};

template <template <typename> class T_struct, typename T_value>
struct value_type<T_struct<T_value>> {
  using type = T_struct<typename value_type<T_value>::type>;
};

template <typename T_value>
struct value_type<std::vector<T_value>> {
  using type = std::vector<typename value_type<T_value>::type>;
};

template <typename T_value, int R, int C>
struct value_type<Eigen::Matrix<T_value, R, C>> {
  using type = Eigen::Matrix<typename value_type<T_value>::type, R, C>;
};

// adjoint_of grabs the adjoint of some data structure

template <typename T>
typename value_type<T>::type adjoint_of(const T x) {
  return typename value_type<T>::type(0);
}

template <typename T>
typename value_type<std::vector<T>>::type adjoint_of(const std::vector<T>& x) {
  const size_t size = x.size();
  typename value_type<std::vector<T>>::type result(size);
  for (size_t i = 0; i != size; ++i) {
    result[i] = adjoint_of(x[i]);
  }
  return result;
}

double adjoint_of(const var& x) { return x.vi_->adj_; }

template <typename T, int R, int C>
inline Eigen::Matrix<typename value_type<T>::type, R, C> adjoint_of(
    const Eigen::Matrix<T, R, C>& M) {
  Eigen::Matrix<typename value_type<T>::type, R, C> Md(M.rows(), M.cols());
  for (int j = 0; j < M.cols(); j++) {
    for (int i = 0; i < M.rows(); i++) {
      Md(i, j) = adjoint_of(M(i, j));
    }
  }
  return Md;
}

template <int R, int C>
inline const Eigen::Matrix<double, R, C>& adjoint_of(
    const Eigen::Matrix<double, R, C>& x) {
  return Eigen::Matrix<double, R, C>::Zero(x.rows(), x.cols());
}

template <int R, int C>
inline const Eigen::Matrix<int, R, C>& adjoint_of(
    const Eigen::Matrix<int, R, C>& x) {
  return Eigen::Matrix<int, R, C>::Zero(x.rows(), x.cols());
}

// as_value is similar to value_of... with the difference that it
// keeps nested container structures which is not the case for
// value_of

template <typename T>
typename value_type<T>::type as_value(const T x) {
  return typename value_type<T>::type(x);
}

template <typename T>
typename value_type<std::vector<T>>::type as_value(const std::vector<T>& x) {
  const size_t size = x.size();
  typename value_type<std::vector<T>>::type result(size);
  for (size_t i = 0; i != size; ++i) {
    result[i] = as_value(x[i]);
  }
  return result;
}

double as_value(const var& x) { return x.vi_->val_; }

template <typename T, int R, int C>
inline Eigen::Matrix<typename value_type<T>::type, R, C> as_value(
    const Eigen::Matrix<T, R, C>& M) {
  Eigen::Matrix<typename value_type<T>::type, R, C> Md(M.rows(), M.cols());
  for (int j = 0; j < M.cols(); j++) {
    for (int i = 0; i < M.rows(); i++) {
      Md(i, j) = as_value(M(i, j));
    }
  }
  return Md;
}

template <int R, int C>
inline const Eigen::Matrix<double, R, C>& as_value(
    const Eigen::Matrix<double, R, C>& x) {
  return x.val();
}

template <int R, int C>
inline const Eigen::Matrix<int, R, C>& as_value(
    const Eigen::Matrix<int, R, C>& x) {
  return x.val();
}

namespace internal {

template <typename T>
inline void add_adjoint(T& sum, const T term) {
  sum += term;
}

template <typename T>
inline void add_adjoint(std::vector<T>& sum, const std::vector<T>& term,
                        int offset = 0) {
  const size_t size = term.size();
  for (size_t i = 0; i != size; ++i) {
    add_adjoint<T>(sum[offset + i], term[i]);
  }
}

template <typename T, int R, int C>
inline void add_adjoint(Eigen::Matrix<T, R, C>& sum,
                        const Eigen::Matrix<T, R, C>& term) {
  sum += term;
}

template <typename Op>
inline void add_adjoints(
    ops_partials_edge<double, std::vector<Op>>& edge,
    const std::vector<typename value_type<Op>::type>& adjoint,
    std::size_t offset = 0) {
  for (std::size_t i = 0; i != adjoint.size(); ++i)
    edge.partials_[offset + i] += adjoint[i];
}

inline vari** register_operands(const var& op, vari** varis) {
  (*varis) = op.vi_;
  return varis + 1;
}

template <int R, int C>
inline vari** register_operands(const Eigen::Matrix<var, R, C>& M_op,
                                vari** varis) {
  const size_t size = M_op.size();
  for (size_t i = 0; i != M_op.cols(); ++i) {
    for (size_t j = 0; j != M_op.rows(); ++j, ++varis) {
      *varis = M_op(i, j).vi_;
    }
  }
  return varis;
}

template <typename T>
inline vari** register_operands(const std::vector<T>& op, vari** varis) {
  const size_t size = op.size();
  for (size_t i = 0; i != size; ++i) {
    varis = register_operands(op[i], varis);
  }
  return varis;
}

inline double* register_partials(const double grad, double* partials) {
  (*partials) = grad;
  return partials + 1;
}

template <int R, int C>
inline double* register_partials(const Eigen::Matrix<double, R, C>& M_grad,
                                 double* partials) {
  const size_t size = M_grad.size();
  std::copy(M_grad.data(), M_grad.data() + size, partials);
  return partials + size;
}

template <typename T>
inline double* register_partials(const std::vector<T>& grad, double* partials) {
  const size_t size = grad.size();
  for (size_t i = 0; i != size; ++i) {
    partials = register_partials(grad[i], partials);
  }
  return partials;
}

template <class ReduceFunction, class M, class T, class Arg1, class Arg2,
          class Arg3, class Arg4>
struct reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2, Arg3, Arg4, var> {
  using vmapped_t = std::vector<M>;
  using arg1_t = std::vector<Arg1>;
  using arg2_t = std::vector<Arg2>;
  using arg3_t = std::vector<Arg3>;
  using arg4_t = std::vector<Arg4>;

  using vmapped_value_t = std::vector<typename value_type<M>::type>;
  using arg1_value_t = std::vector<typename value_type<Arg1>::type>;
  using arg2_value_t = std::vector<typename value_type<Arg2>::type>;
  using arg3_value_t = std::vector<typename value_type<Arg3>::type>;
  using arg4_value_t = std::vector<typename value_type<Arg4>::type>;

  struct recursive_reducer {
    const vmapped_t& vmapped_;
    vmapped_value_t& vmapped_adjoint_;
    const arg1_t& arg1_;
    const arg1_value_t& arg1_value_;
    const arg2_t& arg2_;
    const arg2_value_t& arg2_value_;
    const arg3_t& arg3_;
    const arg3_value_t& arg3_value_;
    const arg4_t& arg4_;
    const arg4_value_t& arg4_value_;

    arg1_value_t arg1_adjoint_;
    arg2_value_t arg2_adjoint_;
    arg3_value_t arg3_adjoint_;
    arg4_value_t arg4_adjoint_;

    double terms_sum_;

    recursive_reducer(const vmapped_t& vmapped,
                      vmapped_value_t& vmapped_adjoint, const T& init,
                      const arg1_t& arg1, const arg1_value_t& arg1_value,
                      const arg2_t& arg2, const arg2_value_t& arg2_value,
                      const arg3_t& arg3, const arg3_value_t& arg3_value,
                      const arg4_t& arg4, const arg4_value_t& arg4_value
                      // ops_partials_t& terms_partials_mapped,
                      )
        : vmapped_(vmapped),
          vmapped_adjoint_(vmapped_adjoint),
          arg1_(arg1),
          arg1_value_(arg1_value),
          arg2_(arg2),
          arg2_value_(arg2_value),
          arg3_(arg3),
          arg3_value_(arg3_value),
          arg4_(arg4),
          arg4_value_(arg4_value),
          arg1_adjoint_(is_constant<arg1_t>::value ? arg1_value_t()
                                                   : adjoint_of(arg1_value)),
          arg2_adjoint_(is_constant<arg2_t>::value ? arg2_value_t()
                                                   : adjoint_of(arg2_value)),
          arg3_adjoint_(is_constant<arg3_t>::value ? arg3_value_t()
                                                   : adjoint_of(arg3_value)),
          arg4_adjoint_(is_constant<arg4_t>::value ? arg4_value_t()
                                                   : adjoint_of(arg4_value)),
          terms_sum_(as_value(init)) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          vmapped_adjoint_(other.vmapped_adjoint_),
          arg1_(other.arg1_),
          arg1_value_(other.arg1_value_),
          arg2_(other.arg2_),
          arg2_value_(other.arg2_value_),
          arg3_(other.arg3_),
          arg3_value_(other.arg3_value_),
          arg4_(other.arg4_),
          arg4_value_(other.arg4_value_),
          arg1_adjoint_(adjoint_of(arg1_value_)),
          arg2_adjoint_(adjoint_of(arg2_value_)),
          arg3_adjoint_(adjoint_of(arg3_value_)),
          arg4_adjoint_(adjoint_of(arg4_value_)),
          terms_sum_(0.0) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      try {
        start_nested();

        // create a deep copy of all var's so that these are not
        // linked to any outer AD tree

        vmapped_t local_sub_slice;

        local_sub_slice.reserve(r.size());
        for (std::size_t i = r.begin(); i != r.end(); ++i)
          local_sub_slice.emplace_back(as_value(vmapped_[i]));

        const arg1_t local_arg1(arg1_value_.begin(), arg1_value_.end());
        const arg2_t local_arg2(arg2_value_.begin(), arg2_value_.end());
        const arg3_t local_arg3(arg3_value_.begin(), arg3_value_.end());
        const arg4_t local_arg4(arg4_value_.begin(), arg4_value_.end());

        T sub_sum_v
            = ReduceFunction()(r.begin(), r.end() - 1, local_sub_slice,
                               local_arg1, local_arg2, local_arg3, local_arg4);

        sub_sum_v.grad();

        terms_sum_ += sub_sum_v.val();

        if (!is_constant_all<M>::value) {
          vmapped_value_t sub_slice_adjoint = adjoint_of(local_sub_slice);
          add_adjoint(vmapped_adjoint_, sub_slice_adjoint, r.begin());
        }

        if (!is_constant<Arg1>::value) {
          add_adjoint(arg1_adjoint_, adjoint_of(local_arg1));
        }
        if (!is_constant<Arg2>::value) {
          add_adjoint(arg2_adjoint_, adjoint_of(local_arg2));
        }
        if (!is_constant<Arg3>::value) {
          add_adjoint(arg3_adjoint_, adjoint_of(local_arg3));
        }
        if (!is_constant<Arg4>::value) {
          add_adjoint(arg4_adjoint_, adjoint_of(local_arg4));
        }
      } catch (const std::exception& e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
    }

    void join(const recursive_reducer& child) {
      terms_sum_ += child.terms_sum_;

      if (!is_constant_all<Arg1>::value) {
        add_adjoint(arg1_adjoint_, child.arg1_adjoint_);
      }
      if (!is_constant_all<Arg2>::value) {
        add_adjoint(arg2_adjoint_, child.arg2_adjoint_);
      }
      if (!is_constant_all<Arg3>::value) {
        add_adjoint(arg3_adjoint_, child.arg3_adjoint_);
      }
      if (!is_constant_all<Arg4>::value) {
        add_adjoint(arg4_adjoint_, child.arg4_adjoint_);
      }
    }
  };

  T operator()(const vmapped_t& vmapped, T init, std::size_t grainsize,
               const arg1_t& arg1, const arg2_t& arg2, const arg3_t& arg3,
               const arg4_t& arg4) const {
    const std::size_t num_jobs = vmapped.size();
    vmapped_value_t vmapped_adjoint = adjoint_of(as_value(vmapped));
    const arg1_value_t arg1_value = as_value(arg1);
    const arg2_value_t arg2_value = as_value(arg2);
    const arg3_value_t arg3_value = as_value(arg3);
    const arg4_value_t arg4_value = as_value(arg4);

    recursive_reducer worker(vmapped, vmapped_adjoint, init, arg1, arg1_value,
                             arg2, arg2_value, arg3, arg3_value, arg4,
                             arg4_value);

    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);

    std::vector<std::size_t> num_terms_arg(5, 0);

    if (!is_constant<vmapped_t>::value)
      num_terms_arg[0] = num_elements(vmapped_adjoint);
    if (!is_constant<Arg1>::value)
      num_terms_arg[1] = num_elements(arg1_value);
    if (!is_constant<Arg2>::value)
      num_terms_arg[2] = num_elements(arg2_value);
    if (!is_constant<Arg3>::value)
      num_terms_arg[3] = num_elements(arg3_value);
    if (!is_constant<Arg4>::value)
      num_terms_arg[4] = num_elements(arg4_value);

    const std::size_t num_terms = sum(num_terms_arg);

    vari** varis
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(num_terms);
    double* partials
        = ChainableStack::instance_->memalloc_.alloc_array<double>(num_terms);

    std::size_t idx = 0;

    if (!is_constant<vmapped_t>::value) {
      register_operands(vmapped, &varis[idx]);
      register_partials(vmapped_adjoint, &partials[idx]);
      idx += num_terms_arg[0];
    }

    if (!is_constant<Arg1>::value) {
      register_operands(arg1, &varis[idx]);
      register_partials(worker.arg1_adjoint_, &partials[idx]);
      idx += num_terms_arg[1];
    }
    if (!is_constant<Arg2>::value) {
      register_operands(arg2, &varis[idx]);
      register_partials(worker.arg3_adjoint_, &partials[idx]);
      idx += num_terms_arg[2];
    }
    if (!is_constant<Arg3>::value) {
      register_operands(arg3, &varis[idx]);
      register_partials(worker.arg3_adjoint_, &partials[idx]);
      idx += num_terms_arg[3];
    }
    if (!is_constant<Arg4>::value) {
      register_operands(arg4, &varis[idx]);
      register_partials(worker.arg4_adjoint_, &partials[idx]);
      idx += num_terms_arg[4];
    }

    return var(new precomputed_gradients_vari(worker.terms_sum_, num_terms,
                                              varis, partials));
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
