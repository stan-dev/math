__kernel void inverse_step1(__global double* ap,__global double* vv,
  int roll_over,int part_size_fixed, int M, int Mpad) {

  int index_0 = get_global_id(0);
  int i = index_0 * part_size_fixed;
  int part_size;
  double faktor;
  if(index_0 < roll_over) {
    i += index_0;
    part_size = part_size_fixed + 1;
  } else {
    i += roll_over;
    part_size = part_size_fixed;
  }

  int offset = index_0 * (part_size_fixed + 1) * (part_size_fixed + 1);

  for(int p = 0; p < part_size; p++) {
    for(int r = 0; r < part_size; r++) {
      if(p == r) {
        vv[offset + p * (part_size_fixed + 1) + r] = 1;
      } else {
        vv[offset + p * (part_size_fixed + 1) + r] = 0;
      }
    }
  }

  for (unsigned int ii = 0; ii < part_size; ii++) {
    if(ii > 0) {
      for (unsigned int j = ii; j < part_size; j++) {
        faktor = ap[(i + ii - 1) * Mpad + (j + i)];
        for (unsigned int k = 0; k < part_size; k++) {
          vv[offset + j * (part_size_fixed + 1) + k]=
            vv[offset + j * (part_size_fixed + 1) + k] - faktor *
            vv[offset + (ii - 1) * (part_size_fixed + 1) + k];
        }
      }
    }

    faktor = ap[(ii + i) * Mpad + ii + i];
    for (unsigned int k = 0; k < part_size; k++) { 
      vv[offset + ii * (part_size_fixed + 1) + k] =
        vv[offset + ii * (part_size_fixed + 1) + k] / faktor;
    }
  }

  for(int p = 0; p < part_size; p++) {
    for(int r = 0; r < part_size; r++) {
      ap[(i + r) * Mpad + (p + i)]=
        vv[offset + p * (part_size_fixed + 1) + r];
    }
  }
}



__kernel void diagonal_mul( __global double * vec1, unsigned int M,
  unsigned int N, unsigned int Npad) { 
  int i = get_global_id(0); 
  int j = get_global_id(1); 
  if(i == j && i < M &&  j < N ){
       vec1[i * Npad + j] *= 0.5;
  }
}
__kernel void copy_lower_tri_upper_tri( __global double * vec1, unsigned int M,
  unsigned int N, unsigned int Npad) { 
  int i = get_global_id(0);
  int j = get_global_id(1);
  if(i < M && j < N && i < j){
       vec1[i * Npad + j] = vec1[j * Npad + i];
  }
}

#define WPT 4
#define RTS  8
#define TS2 32

__kernel void inverse_step2(__global double* ap, __global int* sizes,
  __global double* MM, int repeat, int roll_over,
  int part_size_fixed, int M, int Mpad) { 
 
  int n = get_global_id(2) * 2;
  double sum = 0;
  int part_size1 = 0, part_size2 = 0;
  int offset_i, offset_j;
  
  for(int r = n * repeat; r < (n + 1) * repeat; r++) {
    part_size1 += sizes[r]; 
  } 
  for(int r = (n + 1) * repeat; r < (n + 2) * repeat; r++) {
    part_size2 += sizes[r];
  }
  
  int sizeM = repeat * (part_size_fixed + 1);
  offset_i = (n + 1) * repeat * part_size_fixed;
  offset_j = n * repeat * part_size_fixed;

  if(((n + 1) * repeat) <= roll_over) {
    offset_i += (n + 1) * repeat;
  } else {
    offset_i += roll_over;
  }
  if((n * repeat) <= roll_over) {
    offset_j += n * repeat;
  } else {
    offset_j += roll_over;
  }
  const int row = get_local_id(0);
  const int col = get_local_id(1);
  const int i = TS2 * get_group_id(0) + row;
  const int j = TS2 * get_group_id(1) + col;

  __local double Asub[TS2][TS2];
  __local double Bsub[TS2][TS2];

  double acc[WPT];
  for (int w = 0; w < WPT; w++) {
      acc[w] = 0.0f;
  }

  const int numTiles = (part_size2 + TS2 - 1) / TS2;

  sum = 0;

  for (int t = 0; t < numTiles; t++) {
    for (int w=0; w < WPT; w++) { 
      const int tiledRow = TS2 * t + row; 
      const int tiledCol = TS2 * t + col; 
       
      if(i < part_size2 && (tiledCol + w * RTS) < part_size1) {
        Asub[col + w * RTS][row] = 
          ap[(tiledCol + offset_j + part_size1 + w * RTS) * Mpad + i + offset_i];
      } else {
        Asub[col + w * RTS][row] = 0.0;
      }

      if((j + w * RTS) < part_size1 && tiledRow < part_size2) {
        Bsub[col + w * RTS][row] =
          ap[(j + offset_j + w * RTS) * Mpad + tiledRow + offset_i];
      } else {
        Bsub[col + w * RTS][row] = 0.0;
      }
    }
       
    barrier(CLK_LOCAL_MEM_FENCE); 
   
    for(int k = 0; k < TS2; k++) {         
      for (int w=0; w < WPT; w++) { 
        acc[w] += Asub[k][row] * Bsub[col + w * RTS][k]; 
      }
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  for (int w=0; w<WPT; w++) {
    if(i < part_size2 && (j + w * RTS) < part_size1) {
      MM[(n / 2) * (sizeM) * (sizeM) + i * part_size1 + j + w * RTS] = acc[w];
    }
  }
}

__kernel void inverse_step3(__global double* ap,__global int* sizes,__global double* MM,
  int repeat, int roll_over,int part_size_fixed, int M, int Mpad){ 

  int n = get_global_id(2) * 2;
  double sum = 0;
  int part_size1 = 0, part_size2 = 0;
  int offset_i, offset_j;
  for(int r = n * repeat; r < (n + 1) * repeat; r++) {
    part_size1 += sizes[r];
  }
  for(int r = (n + 1) * repeat; r < (n + 2) * repeat; r++) {
    part_size2 += sizes[r];
  }
  int sizeM = repeat * (part_size_fixed + 1);
  offset_i = (n + 1) * repeat * part_size_fixed;
  offset_j = n * repeat * part_size_fixed;
  if(((n + 1) * repeat) <= roll_over) {
    offset_i += (n + 1) * repeat;
  } else {
    offset_i += roll_over;
  }
  if((n * repeat) <= roll_over) {
    offset_j += n * repeat;
  } else {
    offset_j += roll_over;
  }
  const int row = get_local_id(0);
  const int col = get_local_id(1);
  const int i = TS2 * get_group_id(0) + row;
  const int j = TS2 * get_group_id(1) + col;

  __local double Asub[TS2][TS2];
  __local double Bsub[TS2][TS2];

  double acc[WPT];
  for (int w = 0; w < WPT; w++) {
      acc[w] = 0.0f;
  }

  const int numTiles = (part_size1 + TS2 - 1) / TS2;

  sum = 0;
  for (int t = 0; t < numTiles; t++) {
    for (int w=0; w < WPT; w++) {
      const int tiledRow = TS2 * t + row;
      const int tiledCol = TS2 * t + col;
      if(i < part_size2 && (tiledCol + w * RTS) < part_size1) {
        Asub[col + w * RTS][row] =
          MM[(n / 2) * (sizeM) * (sizeM) + i * part_size1 + tiledCol + w * RTS];
      } else {
        Asub[col + w * RTS][row] = 0.0;
      }
      if((j + w * RTS) < part_size1 && (j + offset_j + w * RTS) < M){
        Bsub[col + w * RTS][row] =
          ap[(j + offset_j + w * RTS) * Mpad + (tiledRow + offset_i - part_size1)];
      }else{
        Bsub[col + w * RTS][row] = 0.0;
      }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int k = 0; k < TS2; k++) {
      for (int w = 0; w < WPT; w++) {
        acc[w] += Asub[k][row] * Bsub[col + w * RTS][k];
      }
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  for (int w = 0; w < WPT; w++) {
   if(i < part_size2 && (j + w * RTS) < part_size1) {
      ap[(j + offset_j + w * RTS) * Mpad + i + offset_i] = -acc[w];
    }
  }   
}  
