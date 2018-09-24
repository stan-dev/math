mutable bool is_finite_ = false;
mutable bool is_cholesky_factor_ = false;
mutable bool is_cholesky_factor_corr_ = false;
mutable bool is_positive_ = false;
mutable bool is_non_negative_ = false;
mutable bool is_corr_ = false;
mutable bool is_cov_ = false;

#ifdef STAN_GPU
mutable cl::Buffer* opencl_buffer_;
#endif
