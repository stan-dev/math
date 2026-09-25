kernel void calculate(__global double* var2_global, double var3, __global int* var4_buffer, __global double* var4_value, __global double* var6_global, int var6_rows, int var6_view, __global double* var7_global, int var7_rows, int var7_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var2 = var2_global[0];
bool var1 = var2 > var3;
bool var4 = var1;
double var6 = 0; if (!((!contains_nonzero(var6_view, LOWER) && j < i) || (!contains_nonzero(var6_view, UPPER) && j > i))) {var6 = var6_global[i + var6_rows * j];}
double var5 = var6 / var2;
var7_global[i + var7_rows * j] = var5;
if(!var4 && atomic_xchg(var4_buffer, 1) == 0){
var4_buffer[1] = i;
var4_buffer[2] = j;
var4_value[0] = var2;
}}
