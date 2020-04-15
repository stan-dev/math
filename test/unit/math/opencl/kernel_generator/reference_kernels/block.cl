kernel void calculate(__global double* var2_global, int var2_rows, int var2_view, int var1_i, int var1_j, __global double* var3_global, int var3_rows, int var3_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && (j + var1_j) < (i + var1_i)) || (!contains_nonzero(var2_view, UPPER) && (j + var1_j) > (i + var1_i)))) {var2 = var2_global[(i + var1_i) + var2_rows * (j + var1_j)];}
double var1 = var2;
var3_global[i + var3_rows * j] = var1;
}
