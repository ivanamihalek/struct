#undef _GLIBCXX_ATOMIC_BUILTINS
#undef _GLIBCXX_USE_INT128

#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/system_error.h>



#ifdef	__cplusplus
extern "C" {
#endif

#include <cuda.h>   
#include "struct.h"

#define TILE_WIDTH 16
#define PITCH 64
#define MEM_SIZE_MAX 134217728L
#define CLOCK_PRECISION  1E9
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/* opt_quat for gpu                             */

/**********************************************************************************/
__device__ int gauss_elim_gpu(float A[4][4],  float solution[4]) {

    float fMaxElem;
    float fAcc;
    float pfVect[4] = {0.0};
    int nDim = 4;
    int i, j, k, m;

    for(k=0; k<(nDim-1); k++){ // base row of matrix
  
	// search of line with max element
	fMaxElem = fabs( A[k][k] );
	m = k;
	for(i=k+1; i<nDim; i++) {
	    if(fMaxElem < fabs(A[i][k]) )  {
		fMaxElem = A[i][k];
		m = i;
	    }
	}
	// permutation of base line (index k) and max element line(index m)
	if(m != k) {
	    for(i=k; i<nDim; i++) {
		fAcc          = A[k][i];
		A[k][i] = A[m][i];
		A[m][i] = fAcc;
	    }
	    fAcc = pfVect[k];
	    pfVect[k] = pfVect[m];
	    pfVect[m] = fAcc;
	}

	if( A[k][k] == 0.) return 1; // needs improvement !!!

	// triangulation of matrix with coefficients
	for(j=(k+1); j<nDim; j++) { // current row of matrix
    
	    fAcc = - A[j][k] / A[k][k];
	    for(i=k; i<nDim; i++) {
		A[j][i] = A[j][i] + fAcc*A[k][i];
	    }
	    pfVect[j] = pfVect[j] + fAcc*pfVect[k]; // free member recalculation
	}
    }
    // in our case the last row should be all zeros
    solution[nDim-1] = 1.0; // arbitrary - we'll renormalize later
    float norm = 1.0;
    for (k=nDim-2; k>=0; k--)  {
	float sol = 0.0;
	for(i= k+1; i<nDim; i++){
	   sol  -= A[k][i]*solution[i];
	}
	sol = sol/A[k][k];
	solution[k] = sol;
	norm += sol*sol;
    }
    norm = sqrt(norm);
    for (k=0; k<nDim; k++) {
	solution[k] /= norm;
    }
    return 0;
}


/**********************************************/
__device__ int swap_gpu (float *a, float *b){
    float c = *a;
    *a = *b;
    *b =  c;
    return 0;
}

/**********************************************/
__device__ int  quadraticSolve_gpu(float *B, float *C,
			       float * root1, float * root2) {
    float discriminant = (*B * *B) - (4 * *C);
      
    //Cannot do imaginary numbers, yet
    if (discriminant < 0) return 1;
      
    float t = -0.5 * ( *B + ((*B < 0) ? -1 : 1) * sqrt(discriminant));
      
    *root1 = t;
    *root2 = *C / t;

    return 0;
}
/**********************************************/
__device__ int  cubicSolve_gpu(float p, float q, float r, 
			   float * root1, float * root2, float * root3) {
    
    float maxSqrt     = sqrt(FLT_MAX);
    
    if (r == 0) {
	//no ant term, so divide by x and the result is a
	//quadratic, but we must include the trivial x = 0 root
	if (quadraticSolve_gpu(&p,&q, root1, root2))
	{
	    *root3 = 0;
	    if (*root1 < *root2) swap_gpu(root1,root2);
	    if (*root2 < 0) 
	    {
		swap_gpu (root2,root3);
		if (*root1 < 0) swap_gpu(root1,root2);
	    }
	    return 3;
	}
	else
	{
	   *root1 = 0;
	    return 1;
	}
    }

    if ((p == 0) && (q == 0))
    {
	//Special case
	//Equation is x^3 == -r
	*root1 = *root2 = *root3 = pow(-r, (float)(1.0/3.0));
	return 3;
    }

    if ((p > maxSqrt) || (p < -maxSqrt))
    {
	//Equation limits to x^3 + p * x^2 == 0
	*root1 = -p;
	return 1;
    }

    if (q > maxSqrt)
    {
	//Special case, if q is large and the root is -r/q,
	//The x^3 term is negligble, and all other terms cancel.
	*root1 = -r / q;
	return 1;
    }

    if (q < -maxSqrt)
    {
	//Special case, equation is x^3 + q x == 0
	*root1 = -sqrt(-q);
	return 1;
    }

    if ((r > maxSqrt) || (r < -maxSqrt))
    {
	//Another special case
	//Equation is x^3 == -r
	*root1 = -pow(r, (float)(1.0/3.0));
	return 1;
    }

    float v = r + (2.0 * p * p / 9.0 - q) * (p / 3.0);

    if ((v > maxSqrt) || (v < -maxSqrt))
    {
	*root1 = -p;
	return 1;
    }

    float uo3 = q / 3.0 - p * p / 9.0;
    float u2o3 = uo3 + uo3;
      
    if ((u2o3 > maxSqrt) || (u2o3 < -maxSqrt))
    {
	if (p==0)
	{
	    if (q > 0)
	    {
		*root1 = -r / q;
		return 1;
	    }
	      
	    if (q < 0)
	    {
		*root1 = -sqrt(-q);
		return 1;
	    }
		
	    *root1 = 0;
	    return 1;
	}

	*root1 = -q/p;
	return 1;
    }

    float uo3sq4 = u2o3 * u2o3;
    if (uo3sq4 > maxSqrt)
    {
	if (p == 0)
	{
	    if (q > 0)
	    {
		*root1 = -r / q;
		return 1;
	    }

	    if (q < 0)
	    {
		*root1 = -sqrt(-q);
		return 1;
	    }

	    *root1 = 0;
	    return 1;
	}

	*root1 = -q / p;
	return 1;
    }

    float j = (uo3sq4 * uo3) + v * v;
  
    if (j > 0)  {//Only one root (but this test can be wrong due to a
	//catastrophic cancellation in j  (i.e. (uo3sq4 * uo3) == v * v)
	float w = sqrt(j);
	if (v < 0)
	    *root1 =         pow((float)(0.5*(w-v)), (float)(1.0/3.0))
 	           - (uo3) * pow((float)(2.0/(w-v)), (float)(1.0/3.0)) - p / 3.0;
	else
	    *root1 = uo3 * pow((float)(2.0/(w+v)), (float)(1.0/3.0))
	         	  -pow((float)(0.5*(w+v)), (float)(1.0/3.0)) - p/3.0;


	return 1;
    }
  
    if (uo3 >= 0) {//Multiple root detected	  
	*root1 = *root2 = *root3 = pow(v, (float)(1.0/3.0)) - p / 3.0;
	return 3;
    }

    float muo3 = - uo3;
    float s;
    if (muo3 > 0)
    {
	s = sqrt(muo3);
	if (p > 0) s = -s;
    }
    else
	s = 0;
      
    float scube = s * muo3;
    if (scube == 0)
    {
	*root1 = - p / 3.0;
	return 1;
    }
      
    float t = - v / (scube + scube);
    float k = acos(t) / 3.0;
    float cosk = cos(k);
    *root1 = (s + s) * cosk - p / 3.0;
      
    float sinsqk = 1.0 - cosk * cosk;
    if (sinsqk < 0) return 1;

    float rt3sink = sqrt(3.0) * sqrt(sinsqk);
    *root2 = s * (-cosk + rt3sink) - p / 3.0;
    *root3 = s * (-cosk - rt3sink) - p / 3.0;

 
    return 3;
}




/**********************************************/
__device__ int quadSolve_gpu(float C, float B, float A, float* root1, float* root2) {
    // Contingency: if A = 0, not a quadratic = linear
    if(A == 0) {
	//If B is zero then we have a NaN
	if(B == 0) return 0;
      
	*root1 = -1.0 * C / B;
	*root2 = *root1;
    }

    float discriminant = (B * B) - (4 * A * C);
      
    //Cannot do imaginary numbers, yet
    if (discriminant < 0) return 0;
      
    float t = -0.5 * ( B + ((B < 0) ? -1 : 1) * sqrt(discriminant));
      
    *root1 = t / A;
    *root2 = C / t;

    return 1;
}
/**************************************************/
__device__ float quarticError_gpu(float a, float b, float c, float d,
		    float roots[4], int  rootCount) {
    float errors[4];
    int root;
    float max_root = FLT_MIN;
    for ( root = 0; root < rootCount; ++ root) {
	
	float value = (((roots[root]+a) * roots[root] + b) * roots[root] + c) * roots[root] + d;

	if (value == 0) { errors[root] = 0; continue; }

	float deriv = ((4 * roots[root] + 3 * a) * roots[root] + 2 * b) * roots[root] + c;
      
	if (deriv != 0) { 
	    errors[root] = abs(value / deriv);
	} else {
	    float secDeriv = (12 * roots[root] + 6 * a) * roots[root] + 2 * b;
	    if (secDeriv != 0)
		errors[root] = sqrt(abs(value / secDeriv));
	    else
	    {
		float thirdDeriv = 24 * roots[root] + 6 * a;
		if (thirdDeriv != 0)
		    errors[root] = pow(abs(value / thirdDeriv), (float)(1.0/3.0));
		else
		    errors[root] = sqrt(sqrt(abs(value)/24));
	    }
	}
	if (max_root < errors[root]) max_root = errors[root];
    }

    return max_root;
}

/**********************************************/
__device__ int ferrariQuarticSolve_gpu(float a, float b, float c, float d,
			float * root1, float * root2, float * root3, float * root4) {
    float rts[4];
    float worst3[3];
    float qrts[4][3]; /* quartic roots for each cubic root */

    if (d == 0.0) {
	*root1 = 0.0;
	return cubicSolve_gpu(a,b,c,root2,root3,root4) + 1;
    }

    int   j;
    int   n4[4];
    float asqinv4;
    float ainv2;
    float d4;
    float yinv2;
    float v1[4],v2[4],v3[4];
    float p,q,r;
    float y;
    float e,f,esq,fsq,ef;
    float g,gg,h,hh;

    ainv2 = a*0.5;
    asqinv4 = ainv2*ainv2;
    d4 = d*4.0 ;

    p = b;
    q = a*c-d4;
    r = (asqinv4 - b)*d4 + c*c;
    int n3 = cubicSolve_gpu (p, q, r, &v3[0],&v3[1],&v3[2]);
    int j3;
    for (j3 = 0; j3 < n3; ++j3)
    {
	y = v3[j3];
	yinv2 = y*0.5;
	esq = asqinv4 - b - y;
	fsq = yinv2*yinv2 - d;
	if ((esq < 0.0) && (fsq < 0.0))
	    n4[j3] = 0;
	else
	{
	    ef = -(0.25*a*y + 0.5*c);
	    if ( ((a > 0.0)&&(y > 0.0)&&(c > 0.0))
		 || ((a > 0.0)&&(y < 0.0)&&(c < 0.0))
		 || ((a < 0.0)&&(y > 0.0)&&(c < 0.0))
		 || ((a < 0.0)&&(y < 0.0)&&(c > 0.0))
		 ||  (a == 0.0)||(y == 0.0)||(c == 0.0))
		/* use ef - */
	    {
		if ((b < 0.0)&&(y < 0.0))
		{
		    e = sqrt(esq);
		    f = ef/e;
		}
		else if (d < 0.0)
		{
		    f = sqrt(fsq);
		    e = ef/f;
		}
		else
		{
		    if (esq > 0.0)
			e = sqrt(esq);
		    else
			e = 0.0;
		    if (fsq > 0.0)
			f = sqrt(fsq);
		    else
			f = 0.0;
		    if (ef < 0.0)
			f = -f;
		}
	    }
	    else
		/* use esq and fsq - */
	    {
		if (esq > 0.0)
		    e = sqrt(esq);
		else
		    e = 0.0;
		if (fsq > 0.0)
		    f = sqrt(fsq);
		else
		    f = 0.0;
		if (ef < 0.0)
		    f = -f;
	    }
	    /* note that e >= 0.0 */
	    g = ainv2 - e;
	    gg = ainv2 + e;
	    if ( ((b > 0.0 && y > 0.0))
		 || ((b < 0.0 && y < 0.0)) )
	    {
		if ((a > 0.0  && e > 0.0)
		    || (a < 0.0  && e < 0.0) )
		    g = (b + y)/gg;
		else
		    if ((a > 0.0  && e < 0.0)
			|| (a < 0.0  && e > 0.0) )
			gg = (b + y)/g;
	    }
	    hh = -yinv2 + f;
	    h = -yinv2 - f;
	    if ( ((f > 0.0 && y < 0.0))
		 || ((f < 0.0 && y > 0.0)) )
		h = d/hh;
	    else
		if ( ((f < 0.0 && y < 0.0))
		     || ((f > 0.0 && y > 0.0)) )
		    hh = d/h;

	    int n1 = quadSolve_gpu(hh, gg, 1.0, &v1[0], &v1[1]);
	    int n2 = quadSolve_gpu(h,   g, 1.0, &v2[0], &v2[1]);
	    n4[j3] = n1*2+n2*2;
	    qrts[0][j3] = v1[0];
	    qrts[1][j3] = v1[1];
	    qrts[n1*2+0][j3] = v2[0];
	    qrts[n1*2+1][j3] = v2[1];
	}
	for (j = 0; j < n4[j3]; ++j)
	    rts[j] = qrts[j][j3];

	worst3[j3] = quarticError_gpu(a, b, c, d, rts, n4[j3]);
    } /* j3 loop */

    j3 = 0;
    if (n3 != 1)
    {
	if ((n4[1] > n4[j3]) ||
	    ((worst3[1] < worst3[j3] ) && (n4[1] == n4[j3]))) j3 = 1;

	if ((n4[2] > n4[j3]) ||
	    ((worst3[2] < worst3[j3] ) && (n4[2] == n4[j3]))) j3 = 2;
    }

    *root1 = qrts[0][j3];
    *root2 = qrts[1][j3];
    *root3 = qrts[2][j3];
    *root4 = qrts[3][j3];

    return (n4[j3]);
}

/*******************************************************************************/
__device__ int construct_ATA_gpu (float ATA[4][4], float a[3], float  b[3]){

    int i,j,k;
    float A[4][4] = {{ 0.0, -b[0], -b[1], -b[2]},
		      {b[0],   0.0, -a[2],  a[1]},
		      {b[1],  a[2],   0.0, -a[0]},
		      {b[2], -a[1],  a[0],   0.0}};
    
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA[i][j] = 0.0;
	    for (k=0; k<4; k++ ) {
		ATA[i][j] += A[k][i]*A[k][j];
	    }
	}
    }
     
    return 0;
}

/*******************************************************************************/
__device__ int add_matrices_gpu  (float matrix1[4][4],float matrix2[4][4],
		   float result[4][4]){
    int i,j;

    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    result[i][j] = matrix1[i][j] + matrix2[i][j];
	}
    }
    return 0;

}

/**********************************************************/
/**********************************************************/
__device__ int opt_quat_gpu(float *x_d, float *x_cm_d, float *y_d, float *y_cm_d,
			    int *set_of_directions_x, int *set_of_directions_y,
			    float q[4], float * rmsd) {

    
    float * x_sub[3];
    float * y_sub[3];
    int  ctr;
    int  i, j;
 
    float ATA     [4][4] = {{0.0}};
    float prev_ATA[4][4] = {{0.0}};
    float ATA_sum [4][4] = {{0.0}};
    float a[3] = {0.0}, b[3] = {0.0};
    
    //int add_matrices  (float matrix1[4][4],float matrix2[4][4],
    //		       float result[4][4]);
    //int construct_ATA (float ATA[4][4], float a[3], float  b[3]);


    if (!3) {
	*rmsd = -1;
	return 1;
    }
    memset ( &(q[0]), 0, 4*sizeof(float) );

    /* find the subset */
    for ( ctr =0; ctr < 3; ctr++ ) {
        i = set_of_directions_x[ctr];
        x_sub[ctr] =  (float*)(((char*)x_d) + (i * PITCH));
        i = set_of_directions_y[ctr];
        y_sub[ctr] =  (float*)(((char*)y_d) + (i * PITCH));
    }  

    /* B = ATA_sum matrix to diagonalize in order to get the quaternion */
    for ( ctr =0; ctr < 3; ctr++ ) {
   	for (i=0; i<3; i++ ) {
	    a[i] = y_sub[ctr][i] + x_sub[ctr][i];
	    b[i] = y_sub[ctr][i] - x_sub[ctr][i];
	}
 	construct_ATA_gpu (ATA, a, b);
	add_matrices_gpu (prev_ATA, ATA, ATA_sum);
	memcpy (prev_ATA[0], ATA_sum[0], 4*4*sizeof(float));
    }
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA_sum[i][j] /= 3;
	}
    }

    /* diagonalize ATA_sum - the eigenvector corresponsing to the
       smallest lambda is the quaternion we are looking for; the
       eigenvalue is the rmsd*/
    float w[4];
    {
	float a,b,c,d,  e,f,g,h,  i,j,k,l,  m,n,o,p; /* synomyms*/
	a =  ATA_sum[0][0]; b =  ATA_sum[0][1]; c =  ATA_sum[0][2]; d =  ATA_sum[0][3]; 
	e =  ATA_sum[1][0]; f =  ATA_sum[1][1]; g =  ATA_sum[1][2]; h =  ATA_sum[1][3]; 
	i =  ATA_sum[2][0]; j =  ATA_sum[2][1]; k =  ATA_sum[2][2]; l =  ATA_sum[2][3]; 
	m =  ATA_sum[3][0]; n =  ATA_sum[3][1]; o =  ATA_sum[3][2]; p =  ATA_sum[3][3]; 

	float B, C, D, E; /* coefficients id the characteristic polynomial - according to sympy */
	/* A = 1.0; will not be used below */
	B = (-a - f - k - p);
	C =  a*f + a*k + a*p - b*e - c*i - d*m + f*k + f*p - g*j - h*n + k*p - l*o;
	D =  - a*f*k - a*f*p + a*g*j + a*h*n - a*k*p + a*l*o + b*e*k + b*e*p - b*g*i
	    - b*h*m - c*e*j + c*f*i + c*i*p - c*l*m - d*e*n + d*f*m - d*i*o + d*k*m
	    - f*k*p + f*l*o + g*j*p - g*l*n - h*j*o + h*k*n;
	E =   a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n - b*e*k*p
	    + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m + c*e*j*p - c*e*l*n
	    - c*f*i*p + c*f*l*m + c*h*i*n - c*h*j*m - d*e*j*o + d*e*k*n + d*f*i*o
	    - d*f*k*m - d*g*i*n + d*g*j*m;

	ferrariQuarticSolve_gpu(B, C, D, E, &w[0], &w[1], &w[2], &w[3]);
    }

    
    float min_w = w[0];
    for (i=1; i<4; i++) {
	if (w[i] > min_w) continue;
	min_w = w[i];
    }
    
    for (i=0; i<4; i++) ATA_sum[i][i] -= min_w;
    if (gauss_elim_gpu (ATA_sum, q)) {
	return 1;
    }
    
    /* rmsd is the lowest eignevalue */
    *rmsd = sqrt(min_w);


    
    return 0;
}



/******************************************************************************/    
/******************************************************************************/    
/******************************************************************************/    
struct greater_rmsd{
    
    __host__ __device__
    bool operator()(Triple x, Triple y) 
    {
        return x.rmsd < y.rmsd;
    }
            
};
    
    
__device__ int  sum_to_zero_gpu (float *a, float *b ) {
    int i;
    float sum,aux;

    sum = 0;
    for (i=0; i<3; i++ ) {
	aux = a[i] + b[i];
	sum += aux*aux;
    }
    return (sum < 0.001);
}


__device__ int  normalized_cross_gpu (float *x, float *y, float * v, float *norm_ptr) {

    /* v is the output */
    float norm = 0;
    float vec[3];
    int i;

    if (sum_to_zero_gpu (x, y) ) return 1;
    
    vec[0] = x[1]*y[2] -  x[2]*y[1];
    norm  += vec[0]*vec[0];
    vec[1] = x[2]*y[0] -  x[0]*y[2]; 
    norm  += vec[1]*vec[1];
    vec[2] = x[0]*y[1] -  x[1]*y[0];
    norm  += vec[2]*vec[2];
    norm   = sqrt (norm);

    for (i=0; i<3; i++ ) {
	v[i] = vec[i] /norm;
    }

    if ( norm_ptr) *norm_ptr = norm;
    

    return 0;
    
}

__device__ int unnorm_dot_gpu (float *x, float *y, float * dot) {


    float cosine = 0;
    int i;
    
    for (i=0; i<3; i++ ) {
	cosine += x[i]*y[i];
    }

    if (cosine > 1.0 ) 
	cosine = 1.0; /* this should be numerical */
   
    *dot = cosine;
   
    return 0;
    
}


__device__ int distance_of_nearest_approach_gpu(float *x_d, float *x_cm_d, float *y_d,
						float *y_cm_d, int *set_of_directions_x,
						int *set_of_directions_y,
        int set_size, float * rmsd_ptr) {

    float cm_vector[3], distance_x, distance_y;
    float aux, rmsd;
    int i, ray_a, ray_b, a, b, prev_next, norm;
    float cross[3];

    //if (set_size <= 1) return 1;
       
    float * row_data_a, * row_data_b;

    
    rmsd = 0.0;
    norm = 0;
    /* the rmsd for the remaining vectors is ... */
    for (a = 0; a < set_size; a++) {

        for (prev_next = -1; prev_next <= +1; prev_next += 2) {
            b = (set_size + a + prev_next) % set_size;

            ray_a = set_of_directions_x[a];
            ray_b = set_of_directions_x[b];
            
            row_data_a = (float*)(((char*)x_cm_d) + (ray_a * PITCH));
            row_data_b = (float*)(((char*)x_cm_d) + (ray_b * PITCH));
    
            
            /* distance of nearest approach of ray b
               to the cm of a, in the set of directions x */
            
            for (i = 0; i < 3; ++i) {
                cm_vector[i] = row_data_b[i] - row_data_a[i];
            }
            
            
            row_data_a = (float*)(((char*)x_d) + (ray_a * PITCH));
            normalized_cross_gpu(cm_vector, row_data_a, cross, &distance_x);

            ray_a = set_of_directions_y[a];
            ray_b = set_of_directions_y[b];
            /* distance of nearest approach of ray b
               to the cm of a, in the set of directions y */
            
            row_data_a = (float*)(((char*)y_cm_d) + (ray_a * PITCH));
            row_data_b = (float*)(((char*)y_cm_d) + (ray_b * PITCH));
    
            for (i = 0; i < 3; ++i) {
                cm_vector[i] = row_data_b[i] - row_data_a[i];
            }
            
            row_data_a = (float*)(((char*)y_d) + (ray_a * PITCH));
            normalized_cross_gpu(cm_vector, row_data_a, cross, &distance_y);

            aux = distance_x - distance_y;
            rmsd += aux*aux;
            norm++;
        }

    }

    rmsd /= norm;
    rmsd = sqrt(rmsd);

    *rmsd_ptr = rmsd;

    return 0;
}
/***************************************************************************/
__device__ int  same_hand_triple_gpu(float *x_d, float *x_cm_d, float *y_d, float *y_cm_d,
				     int *set_of_directions_x, int *set_of_directions_y, int set_size){

    float cm_vector[3], avg_cm[3], cross[3], dx, dy;
    int i, ray_a, ray_b, ray_c;

    // if (set_size != 3) return 0;

    
    /*****************************/
    /*****************************/
    ray_a = set_of_directions_x[0];
    ray_b = set_of_directions_x[1];
    ray_c = set_of_directions_x[2];
    /* I better not use the cross prod here: a small diff
       int he angle makeing them pointing toward each other or
       away from each other changes the direction of cross prod;
       rather, use one vector, and distance between the cm's as the other */
    float * row_data_a, * row_data_b, * row_data_c;

    row_data_a = (float*)(((char*)x_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)x_d) + (ray_b * PITCH));
   
    
    normalized_cross_gpu(row_data_a, row_data_b, cross, NULL);
    
    row_data_a = (float*)(((char*)x_cm_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)x_cm_d) + (ray_b * PITCH));
    row_data_c = (float*)(((char*)x_cm_d) + (ray_c * PITCH)); 
    
    /* note I am making another cm vector here */
    for (i = 0; i < 3; i++) {
        avg_cm[i] = (row_data_b[i] - row_data_a[i])/2;
        cm_vector[i] = row_data_c[i] - avg_cm[i];
        //avg_cm[i] = (x_cm_d[ray_b][i] + x_cm_d[ray_a][i]) / 2;
        //cm_vector[i] = x_cm_d[ray_c][i] - avg_cm[i];
    }
    unnorm_dot_gpu(cm_vector, cross, &dx);

    /*****************************/
    /*****************************/
    ray_a = set_of_directions_y[0];
    ray_b = set_of_directions_y[1];
    ray_c = set_of_directions_y[2];
    
    row_data_a = (float*)(((char*)y_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)y_d) + (ray_b * PITCH));
   
    normalized_cross_gpu(row_data_a, row_data_b, cross, NULL);
    
    row_data_a = (float*)(((char*)y_cm_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)y_cm_d) + (ray_b * PITCH));
    row_data_c = (float*)(((char*)y_cm_d) + (ray_c * PITCH));    
    
    /* note I am making another cm vector here */
    for (i = 0; i < 3; i++) {
        avg_cm[i]    = (row_data_b[i] - row_data_a[i])/2;
        cm_vector[i] = row_data_c[i]- avg_cm[i];
        
        //avg_cm[i] = (y_cm_d[ray_b][i] + y_cm_d[ray_a][i]) / 2;
        //cm_vector[i] = y_cm_d[ray_c][i] - avg_cm[i];
    }
    /*note: unnorm_dot thinks it is getting unit vectors,
      and evrything that is >1 will be "rounded" to 1
      (similarly for -1) - it doesn't do the normalization itself*/
    
    unnorm_dot_gpu(cm_vector, cross, &dy);

    if (dx * dy < 0) return 0;


    return 1; /* this isn't err value - the handedness is the same */
}

/***************************************************************************/

__global__ void find_triplets_gpu (float *x_d, float *x_cm_d, float *y_d, float *y_cm_d,
				   int NX, int NY, int *x_triple_array_d, int *y_triple_array_d, 
				   int cnt_x, int cnt_y, Triple *triple_array_d, float * rmsd_array_d){
    
    int row = blockIdx.x * TILE_WIDTH + threadIdx.x;
    int col = blockIdx.y * TILE_WIDTH + threadIdx.y;
    int i;
    int *row_data;
    float rmsd = -5; 
    float q[4];
    int triple_x[3], triple_y[3];
    // double q_init[4] = {0.0}; // no change

    
    if (row < cnt_x && col < cnt_y) {
	
        triple_array_d[row*cnt_y + col].rmsd = BAD_RMSD + 1;
        rmsd_array_d[row*cnt_y + col] = BAD_RMSD + 1;
        row_data = (int*)(((char*)x_triple_array_d) + (row * PITCH));
        
        for (i = 0; i < 3; ++i) {
            triple_x[i] =  row_data[i];
        }
        
        row_data = (int*)(((char*)y_triple_array_d) + (col * PITCH));
        
        for (i = 0; i < 3; ++i) {
            triple_y[i] =  row_data[i];
        }  
        
        for (i = 0; i <3; ++i) {
            triple_array_d[row*cnt_y + col].triple_x[i] = triple_x[i];
            triple_array_d[row*cnt_y + col].triple_y[i] = triple_y[i];
        }

        if (!same_hand_triple_gpu(x_d, x_cm_d, y_d, y_cm_d, triple_x, triple_y, 3) ) return;

        distance_of_nearest_approach_gpu(x_d, x_cm_d, y_d, y_cm_d, triple_x, triple_y, 3, &rmsd);
        if (rmsd > BAD_RMSD) {
            return; 
        }

	if (opt_quat_gpu (x_d, x_cm_d, y_d, y_cm_d, triple_x, triple_y, q, &rmsd)) return;
	
	// insert values to array of structs
        triple_array_d[row*cnt_y + col].rmsd = rmsd;
        rmsd_array_d[row*cnt_y + col] = rmsd;
        
    }
    
    
}

extern int insert_triple_to_heap_gpu(Representation* X_rep, Representation* Y_rep,
				     int ** x_triple_array, int ** y_triple_array, int x_triple_cnt,
				     int y_triple_cnt, PriorityQueue * heap) {
    

    double **x_db = X_rep->full;
    double **x_cm_db = X_rep->cm;
    double **y_db = Y_rep->full;
    double **y_cm_db = Y_rep->cm;
    int NX = X_rep->N_full;
    int NY = Y_rep->N_full;
    
    float **x = fmatrix(NX, 3);
    float **x_cm = fmatrix(NX, 3);;
    float **y = fmatrix(NY, 3);;
    float **y_cm = fmatrix(NY, 3);;
    
    int i,j;
    
    for(i = 0; i < NX; ++i) {
        for(j = 0; j < 3; ++j) {
            x[i][j] = x_db[i][j];
            x_cm[i][j] = x_cm_db[i][j];
        }
    }
    
    for(i = 0; i < NY; ++i) {
        for(j = 0; j < 3; ++j) {
            y[i][j] = y_db[i][j];
            y_cm[i][j] = y_cm_db[i][j];
        }
    }
    
    float *x_d;
    float *x_cm_d;
    float *y_d;
    float *y_cm_d;
    
    int *x_triple_array_d;
    int *y_triple_array_d;
    
    
    int cnt_x = x_triple_cnt;
    int cnt_y = y_triple_cnt;
    
    struct timespec requestStart, requestEnd;
    clock_gettime(CLOCK_REALTIME, &requestStart);
        
    
    Triple * triple_array = (Triple *) malloc(TOP_RMSD * sizeof(Triple));
    
    
    if (cudaSuccess != cudaMalloc(&x_d, NX * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(x_d,PITCH,x[0],sizeof(float)*3,sizeof(float)*3,NX,cudaMemcpyHostToDevice);   
    
    if (cudaSuccess != cudaMalloc(&x_cm_d, NX * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(x_cm_d,PITCH,x_cm[0],sizeof(float)*3,sizeof(float)*3,NX,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&y_d, NY * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(y_d,PITCH,y[0],sizeof(float)*3,sizeof(float)*3,NY,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&y_cm_d, NY * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(y_cm_d,PITCH,y_cm[0],sizeof(float)*3,sizeof(float)*3,NY,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&x_triple_array_d, cnt_x * 3 * PITCH * sizeof(int))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(x_triple_array_d,PITCH,x_triple_array[0],sizeof(int)*3,sizeof(int)*3,cnt_x,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&y_triple_array_d, cnt_y * 3 * PITCH * sizeof(int))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(y_triple_array_d,PITCH,y_triple_array[0],sizeof(int)*3,sizeof(int)*3,cnt_y,cudaMemcpyHostToDevice);  
  
    
    clock_gettime(CLOCK_REALTIME, &requestEnd);
    
    // call gpu and find triplets
    
    Triple * triple_array_d;   
    float * rmsd_array_d;
    
    size_t size = cnt_x * cnt_y * sizeof(Triple);
    size_t size_ratio = sizeof(Triple)/sizeof(float);

    
    if (size < MEM_SIZE_MAX) {
        if (cudaSuccess != cudaMalloc((void **)&triple_array_d, size)) printf("CUDA allocation error!\n");
        if (cudaSuccess != cudaMalloc((void**)&rmsd_array_d, cnt_x * cnt_y * sizeof(float))) printf("CUDA allocation error!\n");
    } else {
        if (cudaSuccess != cudaMalloc((void **)&triple_array_d, MEM_SIZE_MAX)) printf("CUDA allocation error!\n");
        if (cudaSuccess != cudaMalloc((void**)&rmsd_array_d, MEM_SIZE_MAX / size_ratio)) printf("CUDA allocation error!\n");
    }
 
         
    size_t free_mem, total_mem;
    size_t size_curr = size;
    size_t  y_triple_array_pos = 0;
    
    
    int not_enough_memory = 1;
    
    size_t cnt_y_rest = cnt_y;
    
    
    while (not_enough_memory) {
    
        clock_gettime(CLOCK_REALTIME, &requestStart);
        cudaMemGetInfo(&free_mem, &total_mem);
        
        if (size_curr > MEM_SIZE_MAX){
        
            cnt_y = MEM_SIZE_MAX/(cnt_x * sizeof(Triple));
            cnt_y_rest -= cnt_y;
            
            size = MEM_SIZE_MAX;
            size_curr = cnt_x * cnt_y_rest * sizeof(Triple); 
        } else {
            cnt_y = cnt_y_rest;
            size = size_curr;
            not_enough_memory = 0;
        }
        
        int n_blocks_x = cnt_x/TILE_WIDTH + (cnt_x%TILE_WIDTH == 0 ? 0:1);
        int n_blocks_y = cnt_y/TILE_WIDTH + (cnt_y%TILE_WIDTH == 0 ? 0:1);
        dim3 numBlocks(n_blocks_x, n_blocks_y); 
        dim3 threadsPerBlock(TILE_WIDTH, TILE_WIDTH);
        
        int * y_triple_array_d_curr = (int*)(((char*)y_triple_array_d) + (y_triple_array_pos * PITCH));
        
        find_triplets_gpu <<<numBlocks, threadsPerBlock>>> (x_d, x_cm_d, y_d, y_cm_d, NX, NY, x_triple_array_d, y_triple_array_d_curr, 
        cnt_x, cnt_y, triple_array_d, rmsd_array_d);
 
        cudaThreadSynchronize();
       
        
	Triple * triple_array_output_d;
    
        size_t no_of_out_pairs = TOP_RMSD < (cnt_x * cnt_y)? TOP_RMSD : cnt_x * cnt_y; 
        if (cudaSuccess != cudaMalloc((void **)&triple_array_output_d, no_of_out_pairs*sizeof(Triple))) printf("CUDA allocation error!\n");    
        
        thrust::device_vector<int>  indices(cnt_x * cnt_y); 
        thrust::sequence(indices.begin(),indices.end());
        thrust::device_ptr<Triple> structures(triple_array_d);
        thrust::device_ptr<Triple> structures_out(triple_array_output_d);
        thrust::device_ptr<float>rmsd(rmsd_array_d);
         
        //thrust::sort_by_key(rmsd, rmsd + cnt_x * cnt_y, structures);
        thrust::sort_by_key(rmsd, rmsd + cnt_x * cnt_y, indices.begin());
        

        thrust::device_vector<int>::iterator iter = indices.begin() + no_of_out_pairs;
        
        thrust::gather(indices.begin(), iter, structures, structures_out);

        // thrust::sort(structures, structures+ cnt_x * cnt_y, greater_rmsd());
        
        
        
        //cudaMemcpy(triple_array, triple_array_d, no_of_out_pairs * sizeof(Triple), cudaMemcpyDeviceToHost);
        cudaMemcpy(triple_array, triple_array_output_d, no_of_out_pairs * sizeof(Triple), cudaMemcpyDeviceToHost);
/*
        cudaFree(triple_array_d);
        cudaFree(y_triple_array_d);
        cudaFree(rmsd_array_d);
*/
 
        
        int m;
        for(m = 0; m < no_of_out_pairs; ++m) {
            Insert(triple_array[m], *heap);
        }
    
 
        y_triple_array_pos += cnt_y;
        
        cudaFree(triple_array_output_d);
        
    }
    
    
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(x_cm_d);
    cudaFree(y_cm_d);
    cudaFree(x_triple_array_d);
    cudaFree(y_triple_array_d);
    
    cudaFree(triple_array_d);

    cudaFree(rmsd_array_d);
    
    free(triple_array);
    free_fmatrix(x);
    free_fmatrix(x_cm);
    free_fmatrix(y);
    free_fmatrix(y_cm);
   
    return 0;
}

    

    
#ifdef	__cplusplus
}
#endif
