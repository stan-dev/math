kernel void calculate(__global double* var2_global, int var2_rows, int var2_view, int var3, __global int* var4_buffer, __global double* var4_value){
int i = get_global_id(0);
int j = get_global_id(1);
double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) || (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = var2_global[i + var2_rows * j];}
bool var1 = var2 > var3;
bool var4 = var1;
if(!var4 && atomic_xchg(var4_buffer, 1) == 0){
var4_buffer[1] = i;
var4_buffer[2] = j;
var4_value[0] = var2;
}}
