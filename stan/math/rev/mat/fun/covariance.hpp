#ifndef STAN_MATH_REV_MAT_FUN_COVARIANCE_HPP
#define STAN_MATH_REV_MAT_FUN_COVARIANCE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/err/invalid_argument.hpp>


using namespace Eigen;

namespace stan {
  namespace math {
      
      template <typename T>
      Matrix<T,Dynamic,Dynamic> cov(const Matrix<T,Dynamic,Dynamic>& x){

	// Ensure there are enough observations
	if(x.rows() <= 1)
	  invalid_argument<Matrix<T,Dynamic,Dynamic> >
	    ("cov"
	     , "x"
	     , x
	     , "Too few observations to compute covariance"
	     );

	//Center The matrix
	Matrix<T,Dynamic,Dynamic>
	  centered = x.rowwise() - x.colwise().mean();

	//Create an empty matrix to be viewed
	Matrix<T,Dynamic,Dynamic>
	  cmat(Matrix<T,Dynamic,Dynamic>(x.cols(),x.cols()). //initialize matrix
	       setZero().                                    //zero it
	       template selfadjointView<Lower>().            //operate only on lower tri
	       rankUpdate(centered.adjoint(), 1 / double(x.rows() - 1)) //compute XTX/(n-1)
	       );                                                       //copy up to full size

	return( cmat );
			  
    }

    template <typename T>
    Matrix<T,Dynamic,Dynamic> cov(const std::vector<Matrix<T,Dynamic,1> >& x){
      int nrows = x[0].rows();

      //Ensure there are enough observations
      if(nrows <= 1)
	invalid_argument<unsigned int>
	  ("cov"
	   , "x"
	   , nrows
	   , "Too few observations to compute covariance"
	   );	
      
      Matrix<T,Dynamic,Dynamic> mat(nrows, x.size());

      //Append columns provided they have the right dimensionality
      for(unsigned int i = 0; i < x.size(); ++i){
	if(x[i].rows() != nrows)
	  invalid_argument<unsigned int>
	    ("cov"
	     , "x"
	     , x[i].rows()
	     , "Column vectors not all of same length"
	     );
	
	mat.col(i) = x[i].col(0); 
      }

      //Finish the job with the general matrix covariance
      return( cov(mat) );
    }

    template <typename T>
    Matrix<T,Dynamic,Dynamic> cov(const std::vector<Matrix<T,1,Dynamic> >& x){
      int ncols = x[0].cols();

      //Ensure there are enough observations to compute covariance
      if(x.size() <= 1)
	invalid_argument<unsigned int>
	  ("cov"
	   ,"x"
	   , x.size()
	   , "Too few observations to compute covariance"
	   );
      
      Matrix<T,Dynamic,Dynamic> mat(x.size(), ncols);

      //Append rows provided they have the right dimensionality
      for(unsigned int i = 0; i < x.size(); ++i){
	if(x[i].cols() != ncols)
	  invalid_argument<unsigned int>
	    ("cov"
	     , "x"
	     , x[i].cols()
	     , "Row vectors not all of same length"
	     );
	    
	mat.row(i) = x[i].row(0); 
      }

      //Finish the job with the general matrix covariance	
      return( cov(mat) );
    }

  }
}
#endif
