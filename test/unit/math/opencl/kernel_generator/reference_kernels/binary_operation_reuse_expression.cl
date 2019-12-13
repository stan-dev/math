kernel void calculate(__global double* var1_global, int var1_rows, int var1_view, __global double* var2_global, int var2_rows, int var2_view, __global double* var5_global, int var5_rows, int var5_view, __global double* var9_global, int var9_rows, int var9_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var1 = 0; if (!((!contains_nonzero(var1_view, LOWER) && j < i) || (!contains_nonzero(var1_view, UPPER) && j > i))) {var1 = var1_global[i + var1_rows * j];}
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
double var3 = var1 + var2;
double var4 = var3 * var3;
double var5 = 0; if (!((!contains_nonzero(var5_view, LOWER) && j < i) || (!contains_nonzero(var5_view, UPPER) && j > i))) {var5 = var5_global[i + var5_rows * j];}
double var6 = var4 / var5;
double var7 = var3 / var6;
double var8 = var7 * var6;
var9_global[i + var9_rows * j] = var8;
}