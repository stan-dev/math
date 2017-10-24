R"=====(
__kernel void transpose(
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N ){ 
     a[j*M+i] = b[i*N+j];
    } 
};
      
__kernel void copy(
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N ){ 
     b[i*N+j] = a[i*N+j];
    } 
};  

__kernel void zeros(
          __global double *a,
          unsigned int M,
          unsigned int N, 
          unsigned int part) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N ){ 
     if ( part==0 && j<i ){ 
       a[i*N+j] = 0;
     }else if ( part==1 && j>i ){   
       a[i*N+j] = 0;
     }else if ( part==2 ){   
       a[i*N+j] = 0;
     }
    } 
};  
      
__kernel void identity(
          __global double *a,
          unsigned int M,
          unsigned int N) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N ){ 
       if ( i==j ){ 
         a[i*N+j] = 1.0;
       }else{ 
         a[i*N+j] = 0.0;
       } 
    } 
}; 
      
__kernel void copy_triangular(
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N, 
          unsigned int lower_upper) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N){ 
     if ( !lower_upper && i < M && j <=i ){ 
       a[i*N+j] = b[i*N+j];
     }else if (!lower_upper && i<M && j<N){ 
       a[i*N+j]=0; 
     }else if (lower_upper && i<M && j>=i && j<N){ 
       a[i*N+j]=b[i*N+j];
     }else if (lower_upper && i<M && j<i ){ 
       a[i*N+j]=0;
     } 
    } 
};  

__kernel void copy_triangular_transposed(
          __global double *a,
          unsigned int M,
          unsigned int N, 
          unsigned int lower_to_upper) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N ){ 
     if ( lower_to_upper && j>i){ 
       a[j*N+i]=a[i*N+j];  
     }else if ( !lower_to_upper && j>i){ 
       a[i*N+j]=a[j*N+i];  
     } 
    } 
};
      
__kernel void add(
          __global double *c,
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N ){ 
     c[i*N+j] = a[i*N+j]+b[i*N+j];
    } 
}; 
      
__kernel void subtract(
          __global double *c,
          __global double *a,
          __global double *b,
          unsigned int M,
          unsigned int N) 
{ 
    int i = get_global_id(0); 
    int j = get_global_id(1); 
    if ( i < M && j < N ){ 
     c[i*N+j] = a[i*N+j]-b[i*N+j];
    } 
};
)====="
