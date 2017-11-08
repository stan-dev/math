R"=====(
__kernel void cholesky_block(__global double *b,int offset,int M, int n,__global double *V, __global double *d) {	
	int f=get_local_id(0);
	int arrSize=n;
	for (int i = 0; i < n; i++){
		d[i*n+f]=b[(i+offset)*M+f+offset];
		if (i==f)
			V[i*n+f]=1.0;
		else
			V[i*n+f]=0.0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);	
	
	double sig;
	for (int i=0;i<n-1;i++){
		sig=sqrt(d[i*n+i]);
		if (f==0){
			d[i*n+i]=sig;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (f>=(i+1))
			d[f*n+i]=d[f*n+i]/sig;		
		barrier(CLK_LOCAL_MEM_FENCE);
		for (int k=i+1;k<n;k++){
		  if (f>=(i+1))
		  	d[k*n+f]=d[k*n+f]-d[k*n+i]*d[f*n+i];		  
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	
	barrier(CLK_LOCAL_MEM_FENCE);
	if (f==0){
		d[(n-1)*n+n-1]=sqrt(d[(n-1)*n+n-1]);	
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	for (int p=f+1;p<n;p++){
		d[f*n+p]=0;
	}
	
	
	barrier(CLK_LOCAL_MEM_FENCE);	
	//write to matrix
	for (int i = 0; i < n; i++){
		b[(i+offset)*M+f+offset]=d[i*n+f];
	}
	barrier(CLK_LOCAL_MEM_FENCE);	
	double faktor;
	for (int i = 0; i < n; i++){
		if (i>0){
			barrier(CLK_LOCAL_MEM_FENCE);
			for (int j = i; j < n; j++) {
				faktor=d[j*n+i-1];
				V[j*n+f]=V[j*n+f]-faktor*V[(i-1)*n+f];
				barrier(CLK_LOCAL_MEM_FENCE);
			}
		}
		faktor=d[i*n+i];
		V[i*n+f]=V[i*n+f]/faktor;
		barrier(CLK_LOCAL_MEM_FENCE);	
	}	
}

#define TS3	32
__kernel void cholesky_left_update(__global double* l,__global double* ap, __global double* temp,int offset,int block,int M, int max_threads){
	
	
	int k;
	const int row = get_local_id(0);
    const int col = get_local_id(1);
    const int i = TS3*get_group_id(0) + row;
    const int j = TS3*get_group_id(1) + col;
	
    __local double Asub[TS3][TS3];
    __local double Bsub[TS3][TS3];
	
	double sum=0;
	const int numTiles = (block+TS3-1)/TS3;
	
	for (int t=0; t<numTiles; t++) {
	
		const int tiledRow = TS3*t + row;
		const int tiledCol = TS3*t + col;
	
		if (i<max_threads && tiledCol<block && (+offset+block)<M && (offset+tiledCol)<M){
			Asub[col][row] = ap[(i+offset+block)*M+offset+tiledCol];
		}else{
			Asub[col][row] = 0.0;
		}
		
		if (j<block && tiledRow<block){
			Bsub[row][col] = temp[j*block+tiledRow];
		}else{		
			Bsub[row][col] = 0.0;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		for (k=0;k<TS3;k++){
			sum+=Asub[k][row]*Bsub[k][col];			
		}
		barrier(CLK_LOCAL_MEM_FENCE);		
	}
	if (i<max_threads && j<block && i<M && j<M){
		l[i*block+j]=sum;
	}	
}

__kernel void cholesky_mid_update(__global double* l,__global double* ap,int offset,int block,int M, int max_threads){
	
	int k;
	const int row = get_local_id(0);
   const int col = get_local_id(1);
    const int i = TS3*get_group_id(0) + row;
    const int j = TS3*get_group_id(1) + col;

    __local double Asub[TS3][TS3];
    __local double Bsub[TS3][TS3];
	
	double sum;
	sum=0;
		
	const int numTiles = (block+TS3-1)/TS3;
	for (int t=0; t<numTiles; t++) {
		
		const int tiledRow = TS3*t + row;
		const int tiledCol = TS3*t + col;
		
		if (i<max_threads && tiledCol<max_threads && i<M && tiledCol<M){
			Asub[col][row] = l[i*block+tiledCol];
		}else{
			Asub[col][row] = 0.0;
		}
		
		if (j<max_threads && tiledRow<block && j<M && tiledRow<M){
			Bsub[row][col] = l[j*block+tiledRow ];
		}else{		
			Bsub[row][col] = 0.0;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (i<max_threads && j<max_threads && i<M && j<M){
			for (k=0;k<TS3;k++){
				sum+=Asub[k][row]*Bsub[k][col];
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (i<max_threads && j<max_threads && i<M && j<M){
		ap[(i+offset+block)*M+offset+block+j]-=sum;	
	}	
	if (i<max_threads && j<block){
		ap[(i+offset+block)*M+offset+j]=l[i*block+j];
	}
}

__kernel void cholesky_zero(__global double* A,int M){
		
		
	const int i = get_global_id(0);	
    const int j = get_global_id(1);	
    	
    if (j>i){	
		A[i*M+j]=0.0;	
    }
}
)====="
