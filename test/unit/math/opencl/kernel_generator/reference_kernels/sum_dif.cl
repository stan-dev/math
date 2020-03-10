kernel void calculate(__global double* var2_global, int var2_rows, int var2_view, __global double* var3_global, int var3_rows, int var3_view, __global double* var4_global, int var4_rows, int var4_view, __global double* var6_global, int var6_rows, int var6_view, __global double* var7_global, int var7_rows, int var7_view, __global double* var8_global, int var8_rows, int var8_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
double var3 = 0; if (!((!contains_nonzero(var3_view, LOWER) && j < i) || (!contains_nonzero(var3_view, UPPER) && j > i))) {var3 = var3_global[i + var3_rows * j];}
double var1 = var2 + var3;
var4_global[i + var4_rows * j] = var1;
double var6 = 0; if (!((!contains_nonzero(var6_view, LOWER) && j < i) || (!contains_nonzero(var6_view, UPPER) && j > i))) {var6 = var6_global[i + var6_rows * j];}
double var7 = 0; if (!((!contains_nonzero(var7_view, LOWER) && j < i) || (!contains_nonzero(var7_view, UPPER) && j > i))) {var7 = var7_global[i + var7_rows * j];}
double var5 = var6 - var7;
var8_global[i + var8_rows * j] = var5;
}
