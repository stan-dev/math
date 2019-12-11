kernel void calculate(__global double* var1_global, int var1_rows, int var1_view, __global double* var2_global, int var2_rows, int var2_view, __global double* var4_global, int var4_rows, int var4_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var1 = 0; if (!((!contains_nonzero(var1_view, LOWER) && j < i) || (!contains_nonzero(var1_view, UPPER) && j > i))) {var1 = var1_global[i + var1_rows * j];}
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
double var3 = var1 + var2;
var4_global[i + var4_rows * j] = var3;
}