kernel void calculate(__global double* var2_global, int var2_rows, int var2_view, __global double* var3_global, int var3_rows, int var3_view, const int rows, const int cols){
const int gid_i = get_global_id(0);
const int lid_i = get_local_id(0);
const int lsize_i = get_local_size(0);
const int wg_id_i = get_group_id(0);
const int wg_id_j = get_group_id(1);
const int n_groups_i = get_num_groups(0);
const int blocks_rows = (rows + lsize_i - 1) / lsize_i;
const int blocks_cols = (cols + lsize_i - 1) / lsize_i;
const int i0 = lsize_i * wg_id_i;
const int i = i0 + lid_i;
const int j0 = lsize_i * wg_id_j;
for(int lid_j = 0; lid_j < min(cols - j0, lsize_i); lid_j++){
const int j = j0 + lid_j;
double var1 = 0;
__local double var1_local[LOCAL_SIZE_];
if(i < rows){
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
var1 = var2;
}
var1_local[lid_i] = var1;
barrier(CLK_LOCAL_MEM_FENCE);
for (int step = lsize_i / REDUCTION_STEP_SIZE; step > 0; step /= REDUCTION_STEP_SIZE) {
  if (lid_i < step) {
    for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {
      var1_local[lid_i] = var1_local[lid_i] + var1_local[lid_i + step * i];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
}
if (lid_i == 0) {
var3_global[j * blocks_rows + wg_id_i] = var1_local[0];
}
}
}
