kernel void calculate(__global double* var2_global, int var2_rows, int var2_view, __global double* var3_global, int var3_rows, int var3_view, const int rows, const int cols){
int gid_i = get_global_id(0);
int gid_j = get_global_id(1);
int gsize_i = get_global_size(0);
int gsize_j = get_global_size(1);
for(int j = gid_j; j < cols; j += gsize_j){
for(int i = gid_i; i < rows; i += gsize_i){
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
double var1 = acosh(var2);
var3_global[i + var3_rows * j] = var1;
}
}
}
