kernel void calculate(__global double* var4_global, __global double* var5_global, __global double* var6_global, int var6_rows, int var6_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var4 = var4_global[0];
double var3 = exp((double)var4);
double var5 = var5_global[0];
double var2 = var3 + var5;
var6_global[i + var6_rows * j] = var2;
}
