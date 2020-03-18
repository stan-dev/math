kernel void calculate(__global double* var3_global, int var3_rows, int var3_view, int var4, __global double* var7_global, int var7_rows, int var7_view, int var8, __global double* var9_global, int var9_rows, int var9_view, const int rows, const int cols){
const int gid_i = get_global_id(0);
const int lid_i = get_local_id(0);
const int lsize_i = get_local_size(0);
const int wg_id_i = get_group_id(0);
const int n_groups_i = get_num_groups(0);
const int blocks_rows = (rows + lsize_i - 1) / lsize_i;
const int blocks_cols = (cols + lsize_i - 1) / lsize_i;
for (int idx = wg_id_i; idx < blocks_rows * blocks_cols; idx += n_groups_i){
const int i0 = lsize_i * (idx % blocks_rows);
const int i = i0 + lid_i;
const int j0 = lsize_i * (idx / blocks_rows);
__local double var5_local[LOCAL_SIZE_ * (LOCAL_SIZE_ + 1)];
for(int lid_j = 0; lid_j < min(rows - i0, lsize_i); lid_j++){
const int j = j0 + lid_j;
if(j0 + lid_i < cols){
double var7 = 0; if (!((!contains_nonzero(var7_view, LOWER) && (i - lid_i + lid_j) < (j - lid_j + lid_i)) || (!contains_nonzero(var7_view, UPPER) && (i - lid_i + lid_j) > (j - lid_j + lid_i)))) {var7 = var7_global[(j - lid_j + lid_i) + var7_rows * (i - lid_i + lid_j)];}
double var6 = var7 + var8;
var5_local[lid_i + lid_j * (LOCAL_SIZE_ + 1)] = var6;
}
}
barrier(CLK_LOCAL_MEM_FENCE);
for(int lid_j = 0; lid_j < min(cols - j0, lsize_i); lid_j++){
const int j = j0 + lid_j;
if(i < rows){
double var3 = 0; if (!((!contains_nonzero(var3_view, LOWER) && j < i) || (!contains_nonzero(var3_view, UPPER) && j > i))) {var3 = var3_global[i + var3_rows * j];}
double var2 = var3 + var4;
double var5 = var5_local[lid_j + lid_i * (LOCAL_SIZE_ + 1)];
double var1 = var2 + var5;
var9_global[i + var9_rows * j] = var1;
}
}
}
}
