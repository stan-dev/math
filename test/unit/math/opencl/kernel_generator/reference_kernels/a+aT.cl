kernel void calculate(__global double* var3_global, int var3_rows, int var3_view, int var4, int var8, __global double* var9_global, int var9_rows, int var9_view){
int i = get_global_id(0);
int j = get_global_id(1);
double var3 = 0; if (!((!contains_nonzero(var3_view, LOWER) && j < i) || (!contains_nonzero(var3_view, UPPER) && j > i))) {var3 = var3_global[i + var3_rows * j];}
double var2 = var3 + var4;
double var7 = 0; if (!((!contains_nonzero(var3_view, LOWER) && i < j) || (!contains_nonzero(var3_view, UPPER) && i > j))) {var7 = var3_global[j + var3_rows * i];}
double var6 = var7 + var8;
double var1 = var2 + var6;
var9_global[i + var9_rows * j] = var1;
}
