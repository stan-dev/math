/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen and Alan C. Hindmarsh @ LLNL
 *                Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This header file contains declarations intended for use by
 * generic iterative solvers of Ax = b. The enumeration gives
 * symbolic names for the type  of preconditioning to be used.
 * The function type declarations give the prototypes for the
 * functions to be called within an iterative linear solver, that
 * are responsible for
 *    multiplying A by a given vector v (SUNATimesFn),
 *    setting up a preconditioner P (SUNPSetupFn), and
 *    solving the preconditioner equation Pz = r (SUNPSolveFn).
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_ITERATIVE_H
#define _SUNDIALS_ITERATIVE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*
 * -----------------------------------------------------------------
 * enum : types of preconditioning
 * -----------------------------------------------------------------
 * SUN_PREC_NONE  : The iterative linear solver should not use
 *                  preconditioning.
 *
 * SUN_PREC_LEFT  : The iterative linear solver uses preconditioning
 *                  on the left only.
 *
 * SUN_PREC_RIGHT : The iterative linear solver uses preconditioning
 *                  on the right only.
 *
 * SUN_PREC_BOTH  : The iterative linear solver uses preconditioning
 *                  on both the left and the right.
 * -----------------------------------------------------------------
 */

/* DEPRECATED PREC_NONE: use SUN_PREC_NONE */
/* DEPRECATED PREC_LEFT: use SUN_PREC_LEFT */
/* DEPRECATED PREC_RIGHT: use SUN_PREC_RIGHT */
/* DEPRECATED PREC_BOTH: use SUN_PREC_BOTH */
enum { PREC_NONE, PREC_LEFT, PREC_RIGHT, PREC_BOTH };
enum { SUN_PREC_NONE, SUN_PREC_LEFT, SUN_PREC_RIGHT, SUN_PREC_BOTH };

/*
 * -----------------------------------------------------------------
 * enum : types of Gram-Schmidt routines
 * -----------------------------------------------------------------
 * SUN_MODIFIED_GS  : The iterative solver uses the modified
 *                    Gram-Schmidt routine SUNModifiedGS listed in
 *                    this file.
 *
 * SUN_CLASSICAL_GS : The iterative solver uses the classical
 *                    Gram-Schmidt routine SUNClassicalGS listed in
 *                    this file.
 * -----------------------------------------------------------------
 */

/* DEPRECATED MODIFIED_GS: use SUN_MODIFIED_GS */
/* DEPRECATED CLASSICAL_GS: use SUN_CLASSICAL_GS */
enum { MODIFIED_GS = 1, CLASSICAL_GS = 2 };
enum { SUN_MODIFIED_GS = 1, SUN_CLASSICAL_GS = 2 };

/*
 * -----------------------------------------------------------------
 * Type: SUNATimesFn
 * -----------------------------------------------------------------
 * An SUNATimesFn multiplies Av and stores the result in z. The
 * caller is responsible for allocating memory for the z vector.
 * The parameter A_data is a pointer to any information about A
 * which the function needs in order to do its job. The vector v
 * is unchanged. An SUNATimesFn returns 0 if successful and a
 * non-zero value if unsuccessful.
 * -----------------------------------------------------------------
 */

/* DEPRECATED ATimesFn: use SUNATimesFn */
typedef int (*ATimesFn)(void *A_data, N_Vector v, N_Vector z);
typedef int (*SUNATimesFn)(void *A_data, N_Vector v, N_Vector z);

/*
 * -----------------------------------------------------------------
 * Type: SUNPSetupFn
 * -----------------------------------------------------------------
 * A SUNPSetupFn is an integrator-supplied routine that accesses data
 * stored in the integrator memory structure (P_data), and calls the
 * user-supplied, integrator-specific preconditioner setup routine.
 * -----------------------------------------------------------------
 */

/* DEPRECATED PSetupFn: use SUNPSetupFn */
typedef int (*PSetupFn)(void *P_data);
typedef int (*SUNPSetupFn)(void *P_data);

/*
 * -----------------------------------------------------------------
 * Type: SUNPSolveFn
 * -----------------------------------------------------------------
 * A SUNPSolveFn solves the preconditioner equation Pz = r for the
 * vector z. The caller is responsible for allocating memory for
 * the z vector. The parameter P_data is a pointer to any
 * information about P which the function needs in order to do
 * its job. The parameter lr is input, and indicates whether P
 * is to be taken as the left preconditioner or the right
 * preconditioner: lr = 1 for left and lr = 2 for right.
 * If preconditioning is on one side only, lr can be ignored.
 * If the preconditioner is iterative, then it should strive to
 * solve the preconditioner equation so that
 *     || Pz - r ||_wrms < tol
 * where the weight vector for the WRMS norm may be accessed from
 * the main integrator memory structure.
 * The vector r should not be modified by the SUNPSolveFn.
 * A SUNPSolveFn returns 0 if successful and a non-zero value if
 * unsuccessful.  On a failure, a negative return value indicates
 * an unrecoverable condition, while a positive value indicates
 * a recoverable one, in which the calling routine may reattempt
 * the solution after updating preconditioner data.
 * -----------------------------------------------------------------
 */

/* DEPRECATED PSolveFn: use SUNPSolveFn */
typedef int (*PSolveFn)(void *P_data, N_Vector r, N_Vector z,
                        realtype tol, int lr);
typedef int (*SUNPSolveFn)(void *P_data, N_Vector r, N_Vector z,
                           realtype tol, int lr);

/*
 * -----------------------------------------------------------------
 * Type: SUNQRAddFn
 * -----------------------------------------------------------------
 * A QRAddFn updates a given QR factorization defined by the input
 * parameters:
 *   Q : N_Vector *
 *   R : realtype *
 * with the input vector
 *   f : N_Vector
 *
 * Additional input parameters include:
 *
 *          m : (int) the number of vectors already in the QR factorization
 *
 *       mMax : (int) the maximum number of vectors to be in the QR
 *              factorization (the number of N_Vectors allocated to be in Q)
 *
 * SUNQR_data : (void *) a structure containing any additional inputs
 *              required for the execution of QRAddFn
 *
 * -----------------------------------------------------------------
*/

typedef int (*SUNQRAddFn)(N_Vector *Q, realtype *R, N_Vector f,
                          int m, int mMax, void *QR_data);

/*
 * -----------------------------------------------------------------
 * Function: SUNModifiedGS
 * -----------------------------------------------------------------
 * SUNModifiedGS performs a modified Gram-Schmidt orthogonalization
 * of the N_Vector v[k] against the p unit N_Vectors at
 * v[k-1], v[k-2], ..., v[k-p].
 *
 * v is an array of (k+1) N_Vectors v[i], i=0, 1, ..., k.
 * v[k-1], v[k-2], ..., v[k-p] are assumed to have L2-norm
 * equal to 1.
 *
 * h is the output k by k Hessenberg matrix of inner products.
 * This matrix must be allocated row-wise so that the (i,j)th
 * entry is h[i][j]. The inner products (v[i],v[k]),
 * i=i0, i0+1, ..., k-1, are stored at h[i][k-1]. Here
 * i0=SUNMAX(0,k-p).
 *
 * k is the index of the vector in the v array that needs to be
 * orthogonalized against previous vectors in the v array.
 *
 * p is the number of previous vectors in the v array against
 * which v[k] is to be orthogonalized.
 *
 * new_vk_norm is a pointer to memory allocated by the caller to
 * hold the Euclidean norm of the orthogonalized vector v[k].
 *
 * If (k-p) < 0, then SUNModifiedGS uses p=k. The orthogonalized
 * v[k] is NOT normalized and is stored over the old v[k]. Once
 * the orthogonalization has been performed, the Euclidean norm
 * of v[k] is stored in (*new_vk_norm).
 *
 * SUNModifiedGS returns 0 to indicate success. It cannot fail.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNModifiedGS(N_Vector* v, realtype **h, int k, int p,
                  realtype *new_vk_norm);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNModifiedGS instead")
int ModifiedGS(N_Vector* v, realtype **h, int k, int p,
               realtype *new_vk_norm);

/*
 * -----------------------------------------------------------------
 * Function: SUNClassicalGS
 * -----------------------------------------------------------------
 * SUNClassicalGS performs a classical Gram-Schmidt
 * orthogonalization of the N_Vector v[k] against the p unit
 * N_Vectors at v[k-1], v[k-2], ..., v[k-p]. The parameters v, h,
 * k, p, and new_vk_norm are as described in the documentation
 * for SUNModifiedGS.
 *
 * stemp is a length k+1 array of realtype which can be used as
 * workspace by the SUNClassicalGS routine.
 *
 * vtemp is an N_Vector array of k+1 vectors which can be used as
 * workspace by the SUNClassicalGS routine.
 *
 * SUNClassicalGS returns 0 to indicate success.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNClassicalGS(N_Vector* v, realtype **h, int k, int p,
                   realtype *new_vk_norm, realtype *stemp,
                   N_Vector* vtemp);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNClassicalGS instead")
int ClassicalGS(N_Vector* v, realtype **h, int k, int p,
                realtype *new_vk_norm, realtype *stemp,
                N_Vector* vtemp);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRfact
 * -----------------------------------------------------------------
 * SUNQRfact performs a QR factorization of the Hessenberg matrix H.
 *
 * n is the problem size; the matrix H is (n+1) by n.
 *
 * h is the (n+1) by n Hessenberg matrix H to be factored. It is
 * stored row-wise.
 *
 * q is an array of length 2*n containing the Givens rotations
 * computed by this function. A Givens rotation has the form:
 * | c  -s |
 * | s   c |.
 * The components of the Givens rotations are stored in q as
 * (c, s, c, s, ..., c, s).
 *
 * job is a control flag. If job==0, then a new QR factorization
 * is performed. If job!=0, then it is assumed that the first
 * n-1 columns of h have already been factored and only the last
 * column needs to be updated.
 *
 * SUNQRfact returns 0 if successful. If a zero is encountered on
 * the diagonal of the triangular factor R, then SUNQRfact returns
 * the equation number of the zero entry, where the equations are
 * numbered from 1, not 0. If SUNQRsol is subsequently called in
 * this situation, it will return an error because it could not
 * divide by the zero diagonal entry.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRfact(int n, realtype **h, realtype *q, int job);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNQRFact instead")
int QRfact(int n, realtype **h, realtype *q, int job);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRsol
 * -----------------------------------------------------------------
 * SUNQRsol solves the linear least squares problem
 *
 * min (b - H*x, b - H*x), x in R^n,
 *
 * where H is a Hessenberg matrix, and b is in R^(n+1).
 * It uses the QR factors of H computed by SUNQRfact.
 *
 * n is the problem size; the matrix H is (n+1) by n.
 *
 * h is a matrix (computed by SUNQRfact) containing the upper
 * triangular factor R of the original Hessenberg matrix H.
 *
 * q is an array of length 2*n (computed by SUNQRfact) containing
 * the Givens rotations used to factor H.
 *
 * b is the (n+1)-vector appearing in the least squares problem
 * above.
 *
 * On return, b contains the solution x of the least squares
 * problem, if SUNQRsol was successful.
 *
 * SUNQRsol returns a 0 if successful.  Otherwise, a zero was
 * encountered on the diagonal of the triangular factor R.
 * In this case, SUNQRsol returns the equation number (numbered
 * from 1, not 0) of the zero entry.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRsol(int n, realtype **h, realtype *q, realtype *b);

SUNDIALS_DEPRECATED_EXPORT_MSG("use SUNQRsol instead")
int QRsol(int n, realtype **h, realtype *q, realtype *b);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRAdd_MGS
 * -----------------------------------------------------------------
 * SUNQRAdd_MGS uses Modified Gram Schmidt to update the QR factorization
 * stored in user inputs
 *   - N_Vector *Q
 *   - realtype *R
 * to include the orthonormalized vector input by
 *   - N_Vector df.
 *
 * Additional input parameters include:
 *
 *      m : (int) current number of vectors in QR factorization
 *
 *   mMax : (int) maximum number of vectors that will be in the QR
 *          factorization (the allocated number of N_Vectors in Q)
 *
 * QRdata : (void *) a struct containing any additional temporary
 *          vectors or arrays required for the QRAdd routine
 *
 * On return, Q and R contain the updated Q R factors, if
 * SUNQRAdd_MGS was successful.
 *
 * SUNQRAdd_MGS returns a 0 if successful.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRAdd_MGS(N_Vector *Q, realtype *R, N_Vector df,
                 int m, int mMax, void *QRdata);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRAdd_ICWY
 * -----------------------------------------------------------------
 * SUNQRAdd_ICWY uses the Inverse Compact WY Modified Gram Schmidt
 * method to update the QR factorization stored in user inputs
 *   - N_Vector *Q
 *   - realtype *R
 *   - realtype *T (held within (void *) QRdata)
 * to include the orthonormalized vector input by
 *   - N_Vector df.
 * where the factorization to be updated is of the form
 *   Q * T * R
 *
 * Additional input parameters include:
 *
 *     m :  (int) current number of vectors in QR factorization
 *
 *  mMax : (int) maximum number of vectors that will be in the QR
 *         factorization (the allocated number of N_Vectors in Q)
 *
 * QRdata : (void *) a struct containing any additional temporary
 *          vectors or arrays required for the QRAdd routine
 *
 * QRdata should contain :
 *        N_Vector vtemp, realtype *temp_array (this will be used for T)
 *
 * On return, Q, R, and T contain the updated Q T R factors, if
 * SUNQRAdd_ICWY was successful.
 *
 * SUNQRAdd_ICWY returns a 0 if successful.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRAdd_ICWY(N_Vector *Q, realtype *R, N_Vector df,
                  int m, int mMax, void *QRdata);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRAdd_ICWY_SB
 * -----------------------------------------------------------------
 *  The same function as SUNQRAdd_ICWY but using a single buffer
 *  for global reductions.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRAdd_ICWY_SB(N_Vector *Q, realtype *R, N_Vector df,
                     int m, int mMax, void *QRdata);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRAdd_CGS2
 * -----------------------------------------------------------------
 * SUNQRAdd_CGS2 uses a Classical Gram Schmidt with Reorthogonalization
 * formulation to update the QR factorization stored in user inputs
 *   - N_Vector *Q
 *   - realtype *R
 * to include the orthonormalized vector input by
 *   - N_Vector df.
 *
 * Additional input parameters include:
 *
 *      m : (int) current number of vectors in QR factorization
 *
 *   mMax : (int) maximum number of vectors that will be in the QR
 *          factorization (the allocated number of N_Vectors in Q)
 *
 * QRdata : (void *) a struct containing any additional temporary
 *          vectors or arrays required for the QRAdd routine
 *
 * QRdata should contain :
 *        N_Vector vtemp, N_Vector vtemp2, realtype *temp_array
 *
 * On return, Q and R contain the updated Q R factors, if
 * SUNQRAdd_CGS2 was successful.
 *
 * SUNQRAdd_CGS2 returns a 0 if successful.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRAdd_CGS2(N_Vector *Q, realtype *R, N_Vector df,
                  int m, int mMax, void *QRdata);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRAdd_DCGS2
 * -----------------------------------------------------------------
 * SUNQRAdd_DCGS2 uses a Classical Gram Schmidt with Reorthogonalization
 * formulation that delays reorthogonlization (for the purpose of
 * reducing number of inner products) to update the QR factorization
 * stored in user inputs
 *   - N_Vector *Q
 *   - realtype *R
 * to include the orthonormalized vector input by
 *   - N_Vector df.
 *
 * Additional input parameters include:
 *
 *      m : (int) current number of vectors in QR factorization
 *
 *   mMax : (int) maximum number of vectors that will be in the QR
 *          factorization (the allocated number of N_Vectors in Q)
 *
 * QRdata : (void *) a struct containing any additional temporary
 *          vectors or arrays required for the QRAdd routine
 *
 * QRdata should contain :
 *        N_Vector vtemp, N_Vector vtemp2, realtype *temp_array
 *
 * On return, Q and R contain the updated Q R factors, if
 * SUNQRAdd_DCGS2 was successful.
 *
 * SUNQRAdd_DCGS2 returns a 0 if successful. Otherwise,....
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRAdd_DCGS2(N_Vector *Q, realtype *R, N_Vector df,
                   int m, int mMax, void *QRdata);

/*
 * -----------------------------------------------------------------
 * Function: SUNQRAdd_DCGS2_SB
 * -----------------------------------------------------------------
 *  The same function as SUNQRAdd_DCGS2 but using a single buffer
 *  for global reductions.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT
int SUNQRAdd_DCGS2_SB(N_Vector *Q, realtype *R, N_Vector df,
                      int m, int mMax, void *QRdata);

#ifdef __cplusplus
}
#endif

#endif
