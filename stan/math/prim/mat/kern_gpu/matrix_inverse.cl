R"=====(
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double not supported by OpenCL implementation."
#endif
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision not supported by OpenCL implementation."
#endif

#define A(i,j) A[j*rows+i]
#define V(i,j) V[offset+j*(part_size_fixed+1)+i] 

__kernel void lower_tri_inv_step1(
                __global double* A,
                __global double* V,
                int remainder,
                int part_size_fixed,
                int rows) {
  int indeks = get_global_id(0);
  int i = indeks*part_size_fixed;
  int part_size;
  double faktor;
  if ( indeks < remainder ) {
    i += indeks;
    part_size = part_size_fixed+1;
  } else {
    i += remainder;
    part_size = part_size_fixed;
  }
  int offset = indeks*(part_size_fixed+1)*(part_size_fixed+1);

  for (int p = 0; p < part_size; p++) {
    for (int r = 0; r < part_size; r++) {
      if ( p == r )
        V(p,r) = 1;
      else
        V(p,r) = 0;
    }
  }

  for (unsigned int ii = 0; ii < part_size; ii++) {
    if ( ii > 0 ) {
      for (unsigned int j = ii; j < part_size; j++) {
        faktor = A((j+i),(i+ii-1));
        for (unsigned int k = 0; k < part_size; k++) {
          V(j,k) -= faktor*V((ii-1),k);
        }
      }
    }
    faktor = A((ii+i),(ii+i));
    for (unsigned int k = 0; k < part_size; k++) {
      V(ii,k) /= faktor;
    }
  }
  for (int p = 0; p < part_size; p++) {
    for (int r=0; r < part_size; r++) {
      A((p+i),(i+r)) = V(p,r);
    }
  }
}

#define temp(i,j) temp[(n/2)*(sizeM)*(sizeM)+j*part_size2+i]

#define WPT 4
#define RTS 8
#define TS2 32
__kernel void lower_tri_inv_step2(
                __global double* A,
                __global int* sizes,
                __global double* temp,
                int repeat,
                int remainder,
                int part_size_fixed,
                int rows) {
 int n = get_global_id(2)*2;
 double sum = 0;
 int part_size1 = 0, part_size2 = 0;
 int offset_i, offset_j;
 for (int r = n*repeat; r < (n+1)*repeat; r++)
  part_size1 += sizes[r];
 for (int r = (n+1)*repeat; r < (n+2)*repeat; r++)
  part_size2 += sizes[r];
 int sizeM = repeat*(part_size_fixed+1);
 offset_i = (n+1)*repeat*part_size_fixed;
 offset_j = n*repeat*part_size_fixed;
 if ( ((n+1)*repeat) <= remainder )
  offset_i += (n+1)*repeat;
 else
  offset_i += remainder;
 if ( (n*repeat) <= remainder )
  offset_j += n*repeat;
 else
  offset_j += remainder;
 const int row = get_local_id(0);
 const int col = get_local_id(1);
  const int i = TS2*get_group_id(0) + row;
  const int j = TS2*get_group_id(1) + col;
  __local double Asub[TS2][TS2];
  __local double Bsub[TS2][TS2];
  double acc[WPT];
  for (int w = 0; w < WPT; w++) {
      acc[w] = 0.0;
  }

 const int numTiles = (part_size2+TS2-1)/TS2;
 sum = 0;
 for (int t = 0; t < numTiles; t++) {
  for (int w = 0; w < WPT; w++) {
   const int tiledRow = TS2*t + row;
   const int tiledCol = TS2*t + col;

   if ( i < part_size2 && (tiledCol+w*RTS) < part_size1 &&
      (tiledCol+offset_j+part_size1+w*RTS) < rows  && (i+offset_i) < rows ) {
    Asub[col+w*RTS][row] =
      A((i+offset_i),(tiledCol+offset_j+part_size1+w*RTS));
   } else {
    Asub[col+w*RTS][row] = 0.0;
   }

   if ( (j+w*RTS) < part_size1 && tiledRow < part_size2 &&
        (tiledRow+offset_i) < rows && (j+offset_j+w*RTS) < rows ) {
    Bsub[col+w*RTS][row] = A((tiledRow+offset_i),(j+offset_j+w*RTS));
   } else {
    Bsub[col+w*RTS][row] = 0.0;
   }
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  for (int k = 0; k < TS2; k++) {
   for (int w = 0; w < WPT; w++) {
    acc[w] += Asub[k][row]*Bsub[col+w*RTS][k];
   }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
 }

 for (int w = 0; w < WPT; w++) {
  if ( i < part_size2 && (j+w*RTS) < part_size1 &&
       i < sizeM && (j+w*RTS) < sizeM ) {
   temp(i,(j+w*RTS)) = acc[w];
  }
 }
}

__kernel void lower_tri_inv_step3(
                __global double* A,
                __global int* sizes,
                __global double* temp,
                int repeat,
                int remainder,
                int part_size_fixed,
                int rows) {
 int n = get_global_id(2)*2;
 double sum = 0;
 int part_size1 = 0, part_size2 = 0;
 int offset_i, offset_j;
 for (int r = n*repeat; r < (n+1)*repeat; r++)
  part_size1 += sizes[r];
 for (int r = (n+1)*repeat; r < (n+2)*repeat; r++)
  part_size2 += sizes[r];

 int sizeM = repeat*(part_size_fixed+1);
 offset_i = (n+1)*repeat*part_size_fixed;
 offset_j = n*repeat*part_size_fixed;
 if ( ((n+1)*repeat) <= remainder )
  offset_i += (n+1)*repeat;
 else
  offset_i += remainder;

 if ( (n*repeat) <= remainder )
  offset_j += n*repeat;
 else
  offset_j += remainder;

 const int row = get_local_id(0);
 const int col = get_local_id(1);
  const int i = TS2*get_group_id(0) + row;
  const int j = TS2*get_group_id(1) + col;

  __local double Asub[TS2][TS2];
  __local double Bsub[TS2][TS2];

  double acc[WPT];
  for (int w = 0; w < WPT; w++) {
      acc[w] = 0.0f;
  }

  const int numTiles = (part_size1+TS2-1)/TS2;

 sum = 0;
 for (int t = 0; t < numTiles; t++) {
  for (int w = 0; w < WPT; w++) {
   const int tiledRow = TS2*t + row;
   const int tiledCol = TS2*t + col;
   if ( i < part_size2 && (tiledCol+w*RTS) < part_size1 &&
        i < sizeM && (tiledCol+w*RTS) < sizeM ) {
    Asub[col+w*RTS][row] =
        temp(i,(tiledCol+w*RTS));
   } else {
    Asub[col+w*RTS][row] = 0.0;
   }
   if ( (j+w*RTS) < part_size1 && (j+offset_j+w*RTS) < rows &&
        (tiledRow+offset_i-part_size1) < rows ) {
    Bsub[col+w*RTS][row] =
        A((tiledRow+offset_i-part_size1),(j+offset_j+w*RTS));
   } else {
    Bsub[col+w*RTS][row] = 0.0;
   }
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int k = 0; k < TS2; k++) {
   for (int w = 0; w < WPT; w++) {
    acc[w] += Asub[k][row]*Bsub[col+w*RTS][k];
   }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
 }
 for (int w = 0; w < WPT; w++) {
  if ( i < part_size2 && (j+w*RTS) < part_size1 &&
      (i+offset_i) < rows && (j+offset_j+w*RTS) < rows ) {
   A((i+offset_i),(j+offset_j+w*RTS)) = -acc[w];
  }
 }
}
)====="
