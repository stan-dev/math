kernel void calculate(__global int* var1_global, int var1_rows, int var1_view, __global int* var2_global, int var2_rows, int var2_view, __global double* var3_global, int var3_rows, int var3_view, __global double* var4_global, int var4_rows, int var4_view){
int i = get_global_id(0);
int j = get_global_id(1);
int var1 = 0; if (!((!contains_nonzero(var1_view, LOWER) && j < i) || (!contains_nonzero(var1_view, UPPER) && j > i))) {var1 = var1_global[i + var1_rows * j];}
int var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
double var3 = 0; if (!((!contains_nonzero(var3_view, LOWER) && var2 < var1) || (!contains_nonzero(var3_view, UPPER) && var2 > var1))) {var3 = var3_global[var1 + var3_rows * var2];}
var4_global[i + var4_rows * j] = var3;
}
