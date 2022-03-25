/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * Routines needed to read a sparse matrix in Rutherford-Boeing
 * format. This file is derived from the dreadrb.c in SuperLU_DIST.
 * See original Copyright holder statement below.
 * ----------------------------------------------------------------- */

/*
    Copyright (c) 2003, The Regents of the University of California, through
    Lawrence Berkeley National Laboratory (subject to receipt of any required
    approvals from U.S. Dept. of Energy)

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    (1) Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    (2) Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
    (3) Neither the name of Lawrence Berkeley National Laboratory, U.S. Dept. of
    Energy nor the names of its contributors may be used to endorse or promote
    products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*! @file dreadrb.c
 * \brief Read a matrix stored in Rutherford-Boeing format
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * August 15, 2014
 *
 * </pre>
 *
 * Purpose
 * =======
 *
 * Read a realtype PRECISION matrix stored in Rutherford-Boeing format
 * as described below.
 *
 * Line 1 (A72, A8)
 *      Col. 1 - 72   Title (TITLE)
 *      Col. 73 - 80  Matrix name / identifier (MTRXID)
 *
 * Line 2 (I14, 3(1X, I13))
 *      Col. 1 - 14   Total number of lines excluding header (TOTCRD)
 *      Col. 16 - 28  Number of lines for pointers (PTRCRD)
 *      Col. 30 - 42  Number of lines for row (or variable) indices (INDCRD)
 *      Col. 44 - 56  Number of lines for numerical values (VALCRD)
 *
 * Line 3 (A3, 11X, 4(1X, I13))
 *      Col. 1 - 3    Matrix type (see below) (MXTYPE)
 *      Col. 15 - 28  Compressed Column: Number of rows (NROW)
 *                    Elemental: Largest integer used to index variable (MVAR)
 *      Col. 30 - 42  Compressed Column: Number of columns (NCOL)
 *                    Elemental: Number of element matrices (NELT)
 *      Col. 44 - 56  Compressed Column: Number of entries (NNZERO)
 *                    Elemental: Number of variable indeces (NVARIX)
 *      Col. 58 - 70  Compressed Column: Unused, explicitly zero
 *                    Elemental: Number of elemental matrix entries (NELTVL)
 *
 * Line 4 (2A16, A20)
 *      Col. 1 - 16   Fortran format for pointers (PTRFMT)
 *      Col. 17 - 32  Fortran format for row (or variable) indices (INDFMT)
 *      Col. 33 - 52  Fortran format for numerical values of coefficient matrix
 *                    (VALFMT)
 *                    (blank in the case of matrix patterns)
 *
 * The three character type field on line 3 describes the matrix type.
 * The following table lists the permitted values for each of the three
 * characters. As an example of the type field, RSA denotes that the matrix
 * is real, symmetric, and assembled.
 *
 * First Character:
 *      R Real matrix
 *      C Complex matrix
 *      I integer matrix
 *      P Pattern only (no numerical values supplied)
 *      Q Pattern only (numerical values supplied in associated auxiliary value
 *        file)
 *
 * Second Character:
 *      S Symmetric
 *      U Unsymmetric
 *      H Hermitian
 *      Z Skew symmetric
 *      R Rectangular
 *
 * Third Character:
 *      A Compressed column form
 *      E Elemental form
 *
 * </pre>
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_sparse.h>


#ifdef SUNDIALS_INT64_T
#define IFMT "%ld"
#else
#define IFMT "%d"
#endif

#ifndef DREADRB_DEBUG_LVL
#define DREADRB_DEBUG_LVL 0
#endif

#define DREADRB_ABORT(str) { printf(str); exit(-1); }

/*! \brief Eat up the rest of the current line */
static int DumpLine(FILE *fp)
{
    register int c;
    while ((c = fgetc(fp)) != '\n') ;
    return 0;
}

static int ParseIntFormat(char *buf, sunindextype *num, sunindextype *size)
{
    char *tmp;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp);
    while (*tmp != 'I' && *tmp != 'i') ++tmp;
    ++tmp;
    *size = atoi(tmp);
    return 0;
}

static int ParseFloatFormat(char *buf, sunindextype *num, sunindextype *size)
{
    char *tmp, *period;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
    while (*tmp != 'E' && *tmp != 'e' && *tmp != 'D' && *tmp != 'd'
           && *tmp != 'F' && *tmp != 'f') {
        /* May find kP before nE/nD/nF, like (1P6F13.6). In this case the
           num picked up refers to P, which should be skipped. */
        if (*tmp=='p' || *tmp=='P') {
           ++tmp;
           *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
        } else {
           ++tmp;
        }
    }
    ++tmp;
    period = tmp;
    while (*period != '.' && *period != ')') ++period ;
    *period = '\0';
    *size = atoi(tmp); /*sscanf(tmp, "%2d", size);*/

    return 0;
}

static int ReadVector(FILE *fp, sunindextype n, sunindextype *where, sunindextype perline, sunindextype persize)
{
    register sunindextype i, j, item;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
        fgets(buf, 100, fp);    /* read a line at a time */
        for (j=0; j<perline && i<n; j++) {
            tmp = buf[(j+1)*persize];     /* save the char at that place */
            buf[(j+1)*persize] = 0;       /* null terminate */
            item = atoi(&buf[j*persize]);
            buf[(j+1)*persize] = tmp;     /* recover the char at that place */
            where[i++] = item - 1;
        }
    }

    return 0;
}

static int dReadValues(FILE *fp, sunindextype n, realtype *destination,
        sunindextype perline, sunindextype persize)
{
    register sunindextype i, j, k, s;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
        fgets(buf, 100, fp);    /* read a line at a time */
        for (j=0; j<perline && i<n; j++) {
            tmp = buf[(j+1)*persize];     /* save the char at that place */
            buf[(j+1)*persize] = 0;       /* null terminate */
            s = j*persize;
            for (k = 0; k < persize; ++k) /* No D_ format in C */
                if ( buf[s+k] == 'D' || buf[s+k] == 'd' ) buf[s+k] = 'E';
            destination[i++] = atof(&buf[s]);
            buf[(j+1)*persize] = tmp;     /* recover the char at that place */
        }
    }

    return 0;
}



/*! \brief
 *
 * <pre>
 * On input, nonz/nzval/rowind/colptr represents lower part of a symmetric
 * matrix. On exit, it represents the full matrix with lower and upper parts.
 * </pre>
 */
static void
FormFullA(sunindextype n, sunindextype *nonz, realtype **nzval, sunindextype **rowind, sunindextype **colptr)
{
    register sunindextype i, j, k, col, new_nnz;
    sunindextype *t_rowind, *t_colptr, *al_rowind, *al_colptr, *a_rowind, *a_colptr;
    sunindextype *marker;
    realtype *t_val, *al_val, *a_val;

    al_rowind = *rowind;
    al_colptr = *colptr;
    al_val = *nzval;

    if ( !(marker = (sunindextype *) malloc( (n+1) * sizeof(sunindextype)) ) )
	DREADRB_ABORT("malloc fails for marker[]");
    if ( !(t_colptr = (sunindextype *) malloc( (n+1) * sizeof(sunindextype)) ) )
	DREADRB_ABORT("malloc t_colptr[]");
    if ( !(t_rowind = (sunindextype *) malloc( *nonz * sizeof(sunindextype)) ) )
	DREADRB_ABORT("malloc fails for t_rowind[]");
    if ( !(t_val = (realtype*) malloc( *nonz * sizeof(realtype)) ) )
	DREADRB_ABORT("malloc fails for t_val[]");

    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = al_colptr[j]; i < al_colptr[j+1]; ++i)
	    ++marker[al_rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose matrix A to T */
    for (j = 0; j < n; ++j)
	for (i = al_colptr[j]; i < al_colptr[j+1]; ++i) {
	    col = al_rowind[i];
	    t_rowind[marker[col]] = j;
	    t_val[marker[col]] = al_val[i];
	    ++marker[col];
	}

    new_nnz = *nonz * 2 - n;
    if ( !(a_colptr = (sunindextype *) malloc( (n+1) * sizeof(sunindextype)) ) )
	DREADRB_ABORT("malloc a_colptr[]");
    if ( !(a_rowind = (sunindextype *) malloc( new_nnz * sizeof(sunindextype)) ) )
	DREADRB_ABORT("malloc fails for a_rowind[]");
    if ( !(a_val = (realtype*) malloc( new_nnz * sizeof(realtype)) ) )
	DREADRB_ABORT("malloc fails for a_val[]");

    a_colptr[0] = 0;
    k = 0;
    for (j = 0; j < n; ++j) {
      for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	if ( t_rowind[i] != j ) { /* not diagonal */
	  a_rowind[k] = t_rowind[i];
	  a_val[k] = t_val[i];
	  ++k;
	}
      }

      for (i = al_colptr[j]; i < al_colptr[j+1]; ++i) {
	a_rowind[k] = al_rowind[i];
	a_val[k] = al_val[i];
	++k;
      }

      a_colptr[j+1] = k;
    }

    printf("FormFullA: new_nnz = " IFMT ", k = " IFMT "\n", new_nnz, k);

    free(al_val);
    free(al_rowind);
    free(al_colptr);
    free(marker);
    free(t_val);
    free(t_rowind);
    free(t_colptr);

    *nzval = a_val;
    *rowind = a_rowind;
    *colptr = a_colptr;
    *nonz = new_nnz;
}

void
dreadrb_dist(int iam, FILE *fp, SUNMatrix *Aout, SUNContext sunctx)
{
    register sunindextype i, numer_lines = 0;
    sunindextype tmp, colnum, colsize, rownum, rowsize, valnum, valsize, nrow, ncol, nonz;
    sunindextype *rowind, *colptr;
    realtype *nzval;
    char buf[100], type[4];
    int sym;
    SUNMatrix A;

    /* Line 1 */
    fgets(buf, 100, fp);
    fputs(buf, stdout);

    /* Line 2 */
    for (i=0; i<4; i++) {
        fscanf(fp, "%14c", buf); buf[14] = 0;
        tmp = atoi(buf); /*sscanf(buf, "%d", &tmp);*/
        if (i == 3) numer_lines = tmp;
    }
    DumpLine(fp);

    /* Line 3 */
    fscanf(fp, "%3c", type);
    fscanf(fp, "%11c", buf); /* pad */
    type[3] = 0;
#if (DREADRB_DEBUG_LVL >= 1)
    if ( !iam ) printf("Matrix type %s\n", type);
#endif

    fscanf(fp, "%14c", buf); nrow = atoi(buf);
    fscanf(fp, "%14c", buf); ncol = atoi(buf);
    fscanf(fp, "%14c", buf); nonz = atoi(buf);
    fscanf(fp, "%14c", buf); tmp = atoi(buf);

    if (tmp != 0)
        if ( !iam ) printf("This is not an assembled matrix!\n");
    if (nrow != ncol)
        if ( !iam ) printf("Matrix is not square.\n");
    DumpLine(fp);

    /* Line 4: format statement */
    fscanf(fp, "%16c", buf);
    ParseIntFormat(buf, &colnum, &colsize);
    fscanf(fp, "%16c", buf);
    ParseIntFormat(buf, &rownum, &rowsize);
    fscanf(fp, "%20c", buf);
    ParseFloatFormat(buf, &valnum, &valsize);
    DumpLine(fp);

#if (DREADRB_DEBUG_LVL >= 1)
    if ( !iam ) {
        printf(IFMT " rows, " IFMT " nonzeros\n", *nrow, *nonz);
        printf("colnum " IFMT ", colsize " IFMT "\n", colnum, colsize);
        printf("rownum " IFMT ", rowsize " IFMT "\n", rownum, rowsize);
        printf("valnum " IFMT ", valsize " IFMT "\n", valnum, valsize);
    }
#endif

    A = SUNSparseMatrix(nrow, ncol, nonz, CSC_MAT, sunctx);
    if (A == NULL) DREADRB_ABORT("SUNSparseMatrix returned NULL!\n");

    /* Grab storage for the three arrays ( nzval, rowind, colptr ) */
    nzval = SM_DATA_S(A);
    rowind = SM_INDEXVALS_S(A);
    colptr = SM_INDEXPTRS_S(A);

    nrow = SM_ROWS_S(A);
    ncol = SM_COLUMNS_S(A);
    nonz = SM_NNZ_S(A);

    ReadVector(fp, ncol+1, colptr, colnum, colsize);
    ReadVector(fp, nonz, rowind, rownum, rowsize);
    if ( numer_lines ) {
        dReadValues(fp, nonz, nzval, valnum, valsize);
    }

    sym = (type[1] == 'S' || type[1] == 's');
    if ( sym ) {
	    FormFullA(ncol, &nonz, &nzval, &rowind, &colptr);
    }

    *Aout = A;
}
