#ifndef STAN_MATH_REV_MAT_FUN_CORRELATION_HPP
#define STAN_MATH_REV_MAT_FUN_CORRELATION_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>


using namespace Eigen;

namespace stan {
  namespace math {
      
      template <typename T>
      Matrix<T,Dynamic,Dynamic> corr(const Matrix<T,Dynamic,Dynamic>& x){

	// Ensure there are enough observations
	if(x.rows() <= 1)
	  invalid_argument<Matrix<T,Dynamic,Dynamic> >
	    ("corr"
	     , "x"
	     , x
	     , "Too few observations to compute correlation"
	     );

	//Center The matrix
	Matrix<T,Dynamic,Dynamic>
	  centered = x.rowwise() - x.colwise().mean();

	//Create an empty matrix to be viewed
	SelfAdjointView<Matrix<T,Dynamic,Dynamic>, Lower> cmat =
	  Matrix<T,Dynamic,Dynamic>(x.cols(),x.cols()). //initialize matrix
	  setZero().                                    //zero it
          template selfadjointView<Lower>().            //operate only on lower tri
	  rankUpdate(centered.adjoint(), 1 / double(x.rows() - 1)); //compute XTX/(n-1)
	
	
	std::vector<T> column_variances(x.cols());
	for(int j = 0; j < x.cols(); ++j){
	  Matrix<T,Dynamic,1> curr_col(centered.col(j));
	  column_variances[j] = (curr_col.adjoint() * curr_col)(0) / double(x.rows() - 1);
	}

	for(int j = 0; j < x.cols(); ++j)
	  for(int i = j; i < x.cols(); ++i){
	    if(i == j){
	      cmat(i,j) = column_variances[j];
	    } else {
	      cmat(i,j) /= (sqrt(column_variances[i]) * sqrt(column_variances[j]));
	    }
	  }
	    
	return( Matrix<T,Dynamic,Dynamic>(cmat) );
			  
    }

    template <typename T>
    Matrix<T,Dynamic,Dynamic> corr(const std::vector<Matrix<T,Dynamic,1> >& x){
      int nrows = x[0].rows();

      //Ensure there are enough observations
      if(nrows <= 1)
	invalid_argument<Matrix<T,Dynamic,Dynamic> >
	    ("corr"
	     , "x"
	     , x
	     , "Too few observations to compute correlation"
	     );	
      
      Matrix<T,Dynamic,Dynamic> mat(nrows, x.size());

      //Append columns provided they have the right dimensionality
      for(int i = 0; i < x.size(); ++i){
	if(x[i].rows() != nrows)
	  invalid_argument<Matrix<T,Dynamic,1> >
	    ("corr"
	     , "x"
	     , x
	     , "Column vectors not all of same length"
	     );
	
	mat.col(i) = x[i].col(0); 
      }

      //Finish the job with the general matrix covariance
      return( corr(mat) );
    }

    template <typename T>
    Matrix<T,Dynamic,Dynamic> corr(const std::vector<Matrix<T,1,Dynamic> >& x){
      int ncols = x[0].cols();

      //Ensure there are enough observations to compute covariance
      if(x.size() <= 1)
	invalid_argument<Matrix<T,Dynamic,Dynamic> >
	    ("corr"
	     , "x"
	     , x
	     , "Too few observations to compute correlation"
	     );
      
      Matrix<T,Dynamic,Dynamic> mat(x.size(), ncols);

      //Append rows provided they have the right dimensionality
      for(int i = 0; i < x.size(); ++i){
	if(x[i].cols() != ncols)
	  invalid_argument<Matrix<T,Dynamic,1> >
	    ("corr"
	     , "x"
	     , x
	     , "Column vectors not all of same length"
	     );
	    
	mat.row(i) = x[i].row(0); 
      }

      //Finish the job with the general matrix covariance	
      return( corr(mat) );
    }

  }
}
#endif
