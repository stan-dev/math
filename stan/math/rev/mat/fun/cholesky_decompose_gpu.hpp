#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
// NOTE: for GPU
#include <stan/math/prim/mat/fun/ViennaCL.hpp>
//



    class cholesky_gpu : public vari {
      public:
        int M_;
        vari** variRefA_;
        vari** variRefL_;

        /**
         * Constructor for GPU cholesky function.
         *
         * Stores varis for A.  Instantiates and stores varis for L.
         * Instantiates and stores dummy vari for upper triangular part of var
         * result returned in cholesky_decompose function call
         *
         * variRefL aren't on the chainable autodiff stack, only used for storage
         * and computation. Note that varis for L are constructed externally in
         * cholesky_decompose.
         *
         *
         * @param A matrix
         * @param L_A matrix, cholesky factor of A
         */
        cholesky_gpu(const Eigen::Matrix<var, -1, -1>& A,
                     const Eigen::Matrix<double, -1, -1>& L_A)
          : vari(0.0),
            M_(A.rows()),
            variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                      (A.rows() * (A.rows() + 1) / 2)),
            variRefL_(ChainableStack::memalloc_.alloc_array<vari*>
                      (A.rows() * (A.rows() + 1) / 2)) {
              size_t pos = 0;
              for (size_type j = 0; j < M_; ++j) {
                for (size_type i = j; i < M_; ++i) {
                  variRefA_[pos] = A.coeffRef(i, j).vi_;
                  variRefL_[pos] = new vari(L_A.coeffRef(i, j), false); ++pos;
                }
              }
            }

        /**
         * Reverse mode differentiation algorithm refernce:
         *
         * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
         *
         */
      #ifdef STAN_GPU

        virtual void chain() {

          using Eigen::MatrixXd;
          using Eigen::Lower;
          using Eigen::Block;
          using Eigen::Upper;
          using Eigen::StrictlyUpper;
          using Eigen::StrictlyLower;
          MatrixXd Lbar(M_, M_);
          MatrixXd Lbar1(M_, M_);
          MatrixXd L(M_, M_);
          Lbar.setZero();
          L.setZero();

          int parts;
          if(M_ > 2500){
            parts = 64;
          } else {
            parts = 32;
          }
          //need to setup a smarter way of doing this
          int remainder = M_ % parts;
          int part_size_fixed = M_ / parts;
          int repeat = 1;
          int sizePad;
          
          size_t pos = 0;
          for (size_type j = 0; j < M_; ++j) {
            for (size_type i = j; i < M_; ++i) {
              Lbar.coeffRef(i, j) = variRefL_[pos]->adj_;
              L.coeffRef(i, j) = variRefL_[pos]->val_;
              ++pos;
            }
          }

          viennacl::matrix<double>  vcl_L(M_, M_);
          viennacl::matrix<double>  vcl_temp(M_, M_);
          viennacl::matrix<double>  vcl_Lbar(M_, M_);
          
          viennacl::copy(L, vcl_L);
          viennacl::copy(Lbar, vcl_Lbar);
          vcl_L =  viennacl::trans(vcl_L);
          vcl_Lbar =  viennacl::linalg::prod(vcl_L, vcl_Lbar); 

          my_kernel_mul.global_work_size(0, vcl_Lbar.internal_size1());
          my_kernel_mul.global_work_size(1, vcl_Lbar.internal_size2());		

          my_kernel_mul.local_work_size(0, 32);
          my_kernel_mul.local_work_size(1, 32);

          viennacl::ocl::enqueue(
            my_kernel_mul(
              vcl_Lbar,
              static_cast<cl_uint>(vcl_Lbar.size1()),
              static_cast<cl_uint>(vcl_Lbar.size2()), 
              static_cast<cl_uint>(vcl_Lbar.internal_size2())
            )
          );

          viennacl::matrix<double>  vcl_temp2(M_, M_ * 2);
          std::vector<int> stl_sizes(parts);
          viennacl::vector<int> vcl_sizes(parts);
          
          for(int i = 0; i < parts; i++){
            if(i < remainder)
              stl_sizes[i] = part_size_fixed + 1;
            else
              stl_sizes[i] = part_size_fixed;
          }
          
          viennacl::copy(stl_sizes, vcl_sizes);

          my_kernel_inv1.global_work_size(0, parts);
          my_kernel_inv1.local_work_size(0, 1);

          viennacl::ocl::enqueue(
            my_kernel_inv1(
              vcl_L,vcl_temp,
              static_cast<cl_uint>(remainder),
              static_cast<cl_uint>(part_size_fixed),
              static_cast<cl_uint>(vcl_L.size2()),
              static_cast<cl_uint>(vcl_L.internal_size2())
            )
          );
          

          for(int pp = parts; pp > 1; pp /= 2){

            sizePad = (((part_size_fixed + 1) * repeat + 31) / 32) * 32;
            my_kernel_inv2.global_work_size(0, sizePad);
            my_kernel_inv2.global_work_size(1, sizePad / 4);
            my_kernel_inv2.global_work_size(2, pp / 2);
            my_kernel_inv2.local_work_size(0, 32);
            my_kernel_inv2.local_work_size(1, 32 / 4);
            my_kernel_inv2.local_work_size(2, 1);
            my_kernel_inv3.global_work_size(0, sizePad);
            my_kernel_inv3.global_work_size(1, sizePad / 4);
            my_kernel_inv3.global_work_size(2, pp / 2);
            my_kernel_inv3.local_work_size(0, 32);
            my_kernel_inv3.local_work_size(1, 32/4);
            my_kernel_inv3.local_work_size(2, 1);

            viennacl::ocl::enqueue(
              my_kernel_inv2(
                vcl_L, vcl_sizes, vcl_temp2,
                static_cast<cl_uint>(repeat),
                static_cast<cl_uint>(remainder),
                static_cast<cl_uint>(part_size_fixed),
                static_cast<cl_uint>(vcl_L.size2()),
                static_cast<cl_uint>(vcl_L.internal_size2())
              )
            );
            viennacl::ocl::enqueue(
              my_kernel_inv3(
                vcl_L, vcl_sizes, vcl_temp2,
                static_cast<cl_uint>(repeat),
                static_cast<cl_uint>(remainder),
                static_cast<cl_uint>(part_size_fixed),
                static_cast<cl_uint>(vcl_L.size2()),
                static_cast<cl_uint>(vcl_L.internal_size2())
              )
            );
            repeat *= 2;
          }

          vcl_Lbar =  viennacl::linalg::prod(vcl_L, vcl_Lbar);
          vcl_Lbar =  viennacl::linalg::prod(vcl_L, viennacl::trans(vcl_Lbar));
          viennacl::trans(vcl_Lbar);

          my_kernel_diag.global_work_size(0, vcl_Lbar.internal_size1());
          my_kernel_diag.global_work_size(1, vcl_Lbar.internal_size2());
          my_kernel_diag.local_work_size(0,32);
          my_kernel_diag.local_work_size(1,32);

          viennacl::ocl::enqueue(
            my_kernel_diag(
              vcl_Lbar,
              static_cast<cl_uint>(vcl_Lbar.size1()),
              static_cast<cl_uint>(vcl_Lbar.size2()),
              static_cast<cl_uint>(vcl_Lbar.internal_size2())
            )
          );

          viennacl::copy(vcl_Lbar, Lbar);

          pos = 0;
          for (size_type j = 0; j < M_; ++j)
            for (size_type i = j; i < M_; ++i)
              variRefA_[pos++]->adj_ += Lbar.coeffRef(i, j);

        }
    };
  

#endif  
