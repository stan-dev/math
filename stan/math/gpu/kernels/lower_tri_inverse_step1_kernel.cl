R"(
#ifndef A
#define A(i,j) A[j*rows+i]
#endif

#ifndef V
#define V(i,j) V[offset+j*(part_size_fixed+1)+i] 
#endif

__kernel void lower_tri_inverse_step1(
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
};
)"
