kernel void calculate(__global double* var1_global, int var1_rows, int var1_view, __global double* var3_global, int var3_rows, int var3_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var1 = 0; if (!((!contains_nonzero(var1_view, LOWER) && j < i) || (!contains_nonzero(var1_view, UPPER) && j > i))) {var1 = var1_global[i + var1_rows * j];}
double var2 = acosh(var1);
var3_global[i + var3_rows * j] = var2;
}