kernel void calculate(__global double* var2_global, int var2_rows, int var2_view, __global double* var3_global, int var3_rows, int var3_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
double var1 = acosh((double)var2);
var3_global[i + var3_rows * j] = var1;
}
