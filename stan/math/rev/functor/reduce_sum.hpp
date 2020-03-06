#ifndef STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/fun.hpp>
#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <iostream>
#include <iterator>
#include <tuple>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename ReduceFunction, typename ReturnType, typename Vec,
          typename... Args>
struct reduce_sum_impl<ReduceFunction, require_var_t<ReturnType>, ReturnType,
                       Vec, Args...> {


  struct recursive_reducer {
    size_t num_shared_terms_;
    double* sliced_partials_;
    Vec vmapped_;
    std::ostream* msgs_;
    std::tuple<Args...> args_tuple_;

    double sum_{0.0};
    Eigen::VectorXd args_adjoints_{0};
    template <typename VecT, typename... ArgsT>
    recursive_reducer(size_t num_shared_terms, double* sliced_partials,
                      VecT&& vmapped, std::ostream* msgs,
                      ArgsT&&... args)
        : num_shared_terms_(num_shared_terms),
          sliced_partials_(sliced_partials),
          vmapped_(std::forward<VecT>(vmapped)),
          msgs_(msgs),
          args_tuple_(std::forward<ArgsT>(args)...) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : num_shared_terms_(other.num_shared_terms_),
          sliced_partials_(other.sliced_partials_),
          vmapped_(other.vmapped_),
          msgs_(other.msgs_),
          args_tuple_(other.args_tuple_) {}

    template <typename T, typename = require_arithmetic_t<scalar_type_t<T>>>
    inline auto deep_copy(T arg) {
      return arg;
    }

    inline auto deep_copy(const var& arg) {
      return var(new vari(arg.val(), false));
    }

    template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr>
    inline auto deep_copy(VarVec&& arg) {
      std::vector<var> copy_vec(arg.size());
      for (size_t i = 0; i < arg.size(); ++i) {
        copy_vec[i] = new vari(arg[i].val(), false);
      }
      return copy_vec;
    }

    template <typename VarVec, require_std_vector_st<is_var, VarVec>* = nullptr,
              require_std_vector_vt<is_container, VarVec>* = nullptr>
    inline auto deep_copy(VarVec&& arg) {
      std::vector<value_type_t<VarVec>> copy_vec(arg.size());
      for (size_t i = 0; i < arg.size(); ++i) {
        copy_vec[i] = deep_copy(arg[i]);
      }
      return copy_vec;
    }

    template <typename Mat, require_eigen_vt<is_var, Mat>* = nullptr>
    inline auto deep_copy(Mat&& arg) {
      return arg.unaryExpr([](auto&& x) {
         return var(new vari(x.val(), false));
       }).eval();;
    }

    template <typename... Pargs>
    inline double* accumulate_adjoints(double* dest, const var& x,
                                       Pargs&&... args) {
      *dest += x.adj();
      return accumulate_adjoints(dest + 1, std::forward<Pargs>(args)...);
    }

    template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr, typename... Pargs>
    inline double* accumulate_adjoints(double* dest, VarVec&& x,
                                       Pargs&&... args) {
      for (auto&& x_iter : x) {
        *dest += x_iter.adj();
        ++dest;
      }
      return accumulate_adjoints(dest, std::forward<Pargs>(args)...);
    }

    template <typename VecContainer,
        require_std_vector_st<is_var, VecContainer>* = nullptr,
        require_std_vector_vt<is_container, VecContainer>* = nullptr,
              typename... Pargs>
    inline double* accumulate_adjoints(double* dest, VecContainer&& x,
                                       Pargs&&... args) {
      for (auto&& x_iter : x) {
        dest = accumulate_adjoints(dest, x_iter);
      }
      return accumulate_adjoints(dest, std::forward<Pargs>(args)...);
    }

    template <typename Mat, require_eigen_vt<is_var, Mat>* = nullptr,
              typename... Pargs>
    inline double* accumulate_adjoints(double* dest, const Mat& x,
                                       Pargs&&... args) {
      Eigen::Map<Eigen::MatrixXd>(dest, x.rows(), x.cols()) += x.adj();
      return accumulate_adjoints(dest + x.size(), std::forward<Pargs>(args)...);
    }

    template <typename Arith,
              require_arithmetic_t<scalar_type_t<Arith>>* = nullptr,
              typename... Pargs>
    inline double* accumulate_adjoints(double* dest, Arith&& x,
                                       Pargs&&... args) {
      return accumulate_adjoints(dest, std::forward<Pargs>(args)...);
    }

    inline double* accumulate_adjoints(double* x) { return x; }


    inline void operator()(const tbb::blocked_range<size_t>& r) try {
      if (r.empty())
        return;

      if (args_adjoints_.size() == 0) {
        args_adjoints_ = Eigen::VectorXd::Zero(this->num_shared_terms_);
      }

      start_nested();
      // create a deep copy of all var's so that these are not
      // linked to any outer AD tree
      std::decay_t<Vec> local_sub_slice;
      local_sub_slice.reserve(r.size());
      for (int i = r.begin(); i < r.end(); ++i) {
        local_sub_slice.emplace_back(deep_copy(vmapped_[i]));
      }
      auto args_tuple_local_copy = apply(
          [&](auto&&... args) {
            return std::tuple<decltype(deep_copy(args))...>(
                deep_copy(args)...);
          },
          this->args_tuple_);
      var sub_sum_v = apply(
          [&](auto&&... args) {
            return ReduceFunction()(r.begin(), r.end() - 1, local_sub_slice,
                                    this->msgs_, args...);
          },
          args_tuple_local_copy);
      sub_sum_v.grad();
      sum_ += sub_sum_v.val();
      accumulate_adjoints(this->sliced_partials_ + r.begin(), local_sub_slice);
      apply(
          [&](auto&&... args) {
            return accumulate_adjoints(args_adjoints_.data(),
             std::forward<decltype(args)>(args)...);
          },
          std::move(args_tuple_local_copy));
      recover_memory_nested();
    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }

    void join(const recursive_reducer& rhs) {
      sum_ += rhs.sum_;
      if (args_adjoints_.size() != 0 && rhs.args_adjoints_.size() != 0) {
        args_adjoints_ += rhs.args_adjoints_;
      } else if (args_adjoints_.size() == 0 && rhs.args_adjoints_.size() != 0) {
        args_adjoints_ = rhs.args_adjoints_;
      }
    }
  };

  template <typename... Pargs>
  inline size_t count_var_impl(size_t count, const std::vector<var>& x,
                        Pargs&&... args) const {
    return count_var_impl(count + x.size(), std::forward<Pargs>(args)...);
  }

  template <typename T, require_std_vector_st<is_var, T>* = nullptr,
            require_std_vector_vt<is_container, T>* = nullptr,
            typename... Pargs>
  inline size_t count_var_impl(size_t count, T&& x, Pargs&&... args) const {
    for (size_t i = 0; i < x.size(); i++) {
      count = count_var_impl(count, x[i]);
    }
    return count_var_impl(count, std::forward<Pargs>(args)...);
  }

  template <typename Mat, require_eigen_vt<is_var, Mat>* = nullptr,
            typename... Pargs>
  inline size_t count_var_impl(size_t count, Mat&& x, Pargs&&... args) const {
    return count_var_impl(count + x.size(), std::forward<Pargs>(args)...);
  }

  template <typename... Pargs>
  inline size_t count_var_impl(size_t count, const var& x, Pargs&&... args) const {
    return count_var_impl(count + 1, std::forward<Pargs>(args)...);
  }

  template <typename Arith,
            require_arithmetic_t<scalar_type_t<Arith>>* = nullptr,
            typename... Pargs>
  inline size_t count_var_impl(size_t count, Arith& x, Pargs&&... args) const {
    return count_var_impl(count, std::forward<Pargs>(args)...);
  }

  inline size_t count_var_impl(size_t count) const { return count; }

  /**
   * Count the number of scalars of type T in the input argument list
   *
   * @tparam Pargs Types of input arguments
   * @return Number of scalars of type T in input
   */
  template <typename... Pargs>
  inline size_t count_var(Pargs&&... args) const {
    return count_var_impl(0, std::forward<Pargs>(args)...);
  }

  template <typename... Pargs>
  inline vari** save_varis(vari** dest, const var& x, Pargs&&... args) const {
    *dest = x.vi_;
    return save_varis(dest + 1, std::forward<Pargs>(args)...);
  }

  template <typename VarVec, require_std_vector_vt<is_var, VarVec>* = nullptr,
            typename... Pargs>
  inline vari** save_varis(vari** dest, VarVec&& x, Pargs&&... args) const {
    using write_map = Eigen::Map<Eigen::Matrix<vari*, -1, 1>>;
    using read_map = Eigen::Map<const Eigen::Matrix<var, -1, 1>>;
    write_map(dest, x.size(), 1) = read_map(x.data(), x.size(), 1).vi();
    return save_varis(dest + x.size(), std::forward<Pargs>(args)...);
  }

  template <typename VecContainer,
            require_std_vector_st<is_var, VecContainer>* = nullptr,
            require_std_vector_vt<is_container, VecContainer>* = nullptr,
            typename... Pargs>
  inline vari** save_varis(vari** dest, VecContainer&& x, Pargs&&... args) const {
    for (size_t i = 0; i < x.size(); ++i) {
      dest = save_varis(dest, x[i]);
    }
    return save_varis(dest, std::forward<Pargs>(args)...);
  }

  template <typename Mat, require_eigen_vt<is_var, Mat>* = nullptr,
            typename... Pargs>
  inline vari** save_varis(vari** dest, Mat&& x, Pargs&&... args) const {
    using write_map = Eigen::Map<Eigen::Matrix<vari*,
     std::decay_t<Mat>::RowsAtCompileTime,
     std::decay_t<Mat>::ColsAtCompileTime>>;
    write_map(dest, x.rows(), x.cols()) = x.vi();
    return save_varis(dest + x.size(), std::forward<Pargs>(args)...);
  }

  template <typename R, require_arithmetic_t<scalar_type_t<R>>* = nullptr,
            typename... Pargs>
  inline vari** save_varis(vari** dest, R&& x, Pargs&&... args) const {
    return save_varis(dest, std::forward<Pargs>(args)...);
  }

  inline vari** save_varis(vari** dest) const { return dest; }

  template <typename VecT, typename... ArgsT>
  inline var operator()(VecT&& vmapped, std::size_t grainsize, std::ostream* msgs,
                 ArgsT&&... args) const {
    const std::size_t num_jobs = vmapped.size();

    if (num_jobs == 0)
      return var(0.0);

    const std::size_t num_sliced_terms = count_var(vmapped);
    const std::size_t num_shared_terms = count_var(args...);

    vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
        num_sliced_terms + num_shared_terms);
    double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
        num_sliced_terms + num_shared_terms);

    for (size_t i = 0; i < num_sliced_terms; ++i)
      partials[i] = 0.0;
    recursive_reducer worker(num_shared_terms, partials, vmapped, msgs,
                             args...);

#ifdef STAN_DETERMINISTIC
    tbb::inline_partitioner partitioner;
    tbb::parallel_deterministic_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker,
        partitioner);
#else
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
#endif

    save_varis(varis, std::forward<VecT>(vmapped));
    save_varis(varis + num_sliced_terms, std::forward<ArgsT>(args)...);

    for (size_t i = 0; i < num_shared_terms; ++i) {
      partials[num_sliced_terms + i] = worker.args_adjoints_(i);
    }

    return var(new precomputed_gradients_vari(
        worker.sum_, num_sliced_terms + num_shared_terms, varis, partials));
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
