#include <cblas.h>
#include <omp.h>
#include <sys/time.h>
#include <math.h>

#define n 100
#define iter 100


// timing helper
float timedifference_msec(struct timeval t0, struct timeval t1)
{
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}

// function to add vectors
void vector_add(double* x,
         double a,
         double* y,
         int p){
    
    cblas_daxpy(p,a,x,1,y,1);
    
    }

// function to dot product
double dot(double* a,
         double* b,
         int     p){
    
   return cblas_ddot(p, a, 1, b, 1);
    
}

// copy function wrapper
void copy(double* x,
          double* y,
          int p){
    
    cblas_dcopy(p,x,1,y,1);
    
}

// function for componentwise product
void cprod(double* v, double* w, double* s){
    
    for(int i=0; i<n; i++){
        s[i] = v[i]*w[i];
    }

}


// function to sum entires of vector except entry i

double sum_i(double* v, int i){
    
    double sum = 0.0;
    for(int j=0; j<i; j++){
        sum+= v[j];
    }

    for(int j = (n-1); j>i; j--){
        sum+= v[j];
    }
    return sum;
}

void jacobi(double* A, double* b, double* x0, double* it, double* t, double* error, double* v){
	
    struct timeval t0;
    struct timeval t1;
    float elapsed;

    double s[n];
    double ss[n];
    double x1[n];
    
    copy(x0, x1, n);
    error[0] = 1;
    it[0] = 1;
    for(int k = 1; k < iter; k++){
	it[k] = k+1;	

	gettimeofday(&t0, 0);		
	for(int i = 0; i < n; i++){
	    cprod(A+i*n, x0, s);
            x1[i] = (b[i] - sum_i(s, i))/A[i*n + i];
        }
	copy(x1, x0, n);
	copy(v,ss, n);
	vector_add(x1, -1, ss, n); 
	error[k] = sqrt(dot(ss, ss, n))/sqrt(dot(v,v,n));	
	gettimeofday(&t1, 0);
	elapsed = timedifference_msec(t0, t1);
	t[k] = t[k-1] + elapsed;
    }
    copy(x1, x0, n);
}

void jacobi_parallel(double* A, double* b, double* x0, double* it, double* t, double* error, double* v){
   
    struct timeval t0;
    struct timeval t1;
    float elapsed;

    double s[n];
    double ss[n];
    double x1[n];

    copy(x0, x1, n); 
    error[0] = 1;
    it[0] = 1;
    for(int k = 1; k < iter; k++){
	it[k] = k+1;
	gettimeofday(&t0, 0);

	#pragma omp parallel num_threads(4)
	{
	double s[n]; 
	#pragma omp for
	for(int i = 0; i < n; i++){
	    cprod(A+i*n, x0, s);
            x1[i] = (b[i] - sum_i(s, i))/A[i*n + i];
        }
	}
	copy(x1, x0, n);
	copy(v,ss, n);
	vector_add(x1, -1, ss, n); 
	error[k] = sqrt(dot(ss, ss, n))/sqrt(dot(v,v,n));	
	gettimeofday(&t1, 0);
	elapsed = timedifference_msec(t0, t1);
	t[k] = t[k-1] + elapsed;

    }
    copy(x1, x0, n);
}


void gauss_seidel(double* A, double* b, double* x0, double* it, double* t, double* error, double* v){

    struct timeval t0;
    struct timeval t1;
    float elapsed;

    double s[n];
    double ss[n];

    error[0] = 1;
    it[0] = 1;
    for(int k = 1; k < iter; k++){
	it[k] = k+1;
	gettimeofday(&t0, 0);

	for(int i = 0; i < n; i++){
	    cprod(A+i*n, x0, s);
            x0[i] = (b[i] - sum_i(s, i))/A[i*n + i];
        }
	copy(v,ss, n);
	vector_add(x0, -1, ss, n); 
	error[k] = sqrt(dot(ss, ss, n))/sqrt(dot(v,v,n));	
	gettimeofday(&t1, 0);
	elapsed = timedifference_msec(t0, t1);
	t[k] = t[k-1] + elapsed;

    }
}


