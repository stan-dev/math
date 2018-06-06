R"(
__kernel void lower_tri_inverse_step3(
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
};
)"
