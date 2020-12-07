# ARKODE + XBraid C++ Examples

Examples using ARKODE and XBraid for parallel-in-time integration

## Examples

* `ark_heat2D_xbraid` : 2-D heat equation with SUNDIALS PCG linear solver and
diagonal preconditioner. Serial in space and parallel in time.

* `ark_heat2D_p_xbraid` : 2-D heat equation with SUNDIALS PCG linear solver
and diagonal preconditioner. Parallel in space and time.

* `ark_heat2D_hypre_pfmg_xbraid` : 2-D heat equation with SUNDIALS PCG linear
solver and user-supplied hypre PFMG preconditioner. Parallel in space and time.

## Sample results

The example output files were produced by running:
```
mpiexec -n 4 ./ark_heat2D_xbraid
mpiexec -n 8 ./ark_heat2D_p_xbraid --np 2 2 2
mpiexec -n 8 ./ark_heat2D_hypre_pfmg_xbraid --np 2 2 2
```

The following CMake command was used to configure SUNDIALS:
```
cmake \
 -DCMAKE_INSTALL_PREFIX:PATH="$PWD/../install_opt" \
 -DEXAMPLES_INSTALL_PATH:PATH="$PWD/../examples_opt" \
 -DCMAKE_C_COMPILER=/usr/local/mpich-3.2.1/gnu/bin/mpicc \
 -DMPI_C_COMPILER=/usr/local/mpich-3.2.1/gnu/bin/mpicc \
 -DCMAKE_Fortran_COMPILER=/usr/local/mpich-3.2.1/gnu/bin/mpif90 \
 -DMPI_Fortran_COMPILER=/usr/local/mpich-3.2.1/gnu/bin/mpif90 \
 -DCMAKE_CXX_COMPILER=/usr/local/mpich-3.2.1/gnu/bin/mpicxx \
 -DMPI_CXX_COMPILER=/usr/local/mpich-3.2.1/gnu/bin/mpicxx \
 -DMPIEXEC_EXECUTABLE=/usr/local/mpich-3.2.1/gnu/bin/mpiexec \
 -DCMAKE_C_FLAGS:STRING="-O0 -g -fPIC" \
 -DCMAKE_Fortran_FLAGS:STRING="-O0 -g -fPIC" \
 -DCMAKE_CXX_FLAGS:STRING="-O0 -g -fPIC" \
 -DOpenMP_C_FLAGS:STRING="-fopenmp" \
 -DOpenMP_CXX_FLAGS:STRING="-fopenmp" \
 -DSUNDIALS_PRECISION:STRING="double" \
 -DSUNDIALS_INDEX_TYPE:STRING="int64_t" \
 -DEXAMPLES_ENABLE_C:BOOL="1" \
 -DEXAMPLES_ENABLE_CXX:BOOL="1" \
 -DEXAMPLES_ENABLE_F90:BOOL="1" \
 -DEXAMPLES_INSTALL:BOOL="1" \
 -DFCMIX_ENABLE:BOOL="1" \
 -DMPI_ENABLE:BOOL="1" \
 -DOPENMP_ENABLE:BOOL="1" \
 -DPTHREAD_ENABLE:BOOL="1" \
 -DBUILD_SHARED_LIBS:BOOL="1" \
 -DHYPRE_ENABLE:BOOL="1" \
 -DHYPRE_INCLUDE_DIR:STRING="$HOME/workspace/hypre-2.14.0/gnu_int64/include" \
 -DHYPRE_LIBRARY_DIR:FILEPATH="$HOME/workspace/hypre-2.14.0/gnu_int64/lib" \
 -DXBRAID_DIR:FILEPATH="$HOME/workspace/xbraid" \
 ../sundials
```

Test system information:
System Architecture: x86_64
Processor Type: Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz
Operating System: Ubuntu 18.04
C/Fortran Compilers: gcc/gfortran v7.3.0
MPI: MPICH v3.2.1
