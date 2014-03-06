# include <struct.h>
# include <float.h>

/**********************************************************************************/
int gauss_elim(double A[4][4],  double solution[4]) {

    double fMaxElem;
    double fAcc;
    double pfVect[4] = {0.0};
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
    double norm = 1.0;
    for (k=nDim-2; k>=0; k--)  {
	double sol = 0.0;
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
int swap (float *a, float *b){
    float c = *a;
    *a = *b;
    *b =  c;
    return 0;
}

/**********************************************/
int  quadraticSolve(float *B, float *C,
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
int  cubicSolve(float p, float q, float r, 
	       float * root1, float * root2, float * root3) {
    
    float maxSqrt     = sqrt(FLT_MAX);
    
    if (r == 0) {
	//no ant term, so divide by x and the result is a
	//quadratic, but we must include the trivial x = 0 root
	if (quadraticSolve(&p,&q, root1, root2))
	{
	    *root3 = 0;
	    if (*root1 < *root2) swap(root1,root2);
	    if (*root2 < 0) 
	    {
		swap(root2,root3);
		if (*root1 < 0) swap(root1,root2);
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
	*root1 = *root2 = *root3 = pow(-r, 1.0 / 3.0);
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
	*root1 = -pow(r, 1.0 / 3.0);
	return 1;
    }

    double v = r + (2.0 * p * p / 9.0 - q) * (p / 3.0);

    if ((v > maxSqrt) || (v < -maxSqrt))
    {
	*root1 = -p;
	return 1;
    }

    double uo3 = q / 3.0 - p * p / 9.0;
    double u2o3 = uo3 + uo3;
      
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

    double uo3sq4 = u2o3 * u2o3;
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

    double j = (uo3sq4 * uo3) + v * v;
  
    if (j > 0)  {//Only one root (but this test can be wrong due to a
	//catastrophic cancellation in j  (i.e. (uo3sq4 * uo3) == v * v)
	double w = sqrt(j);
	if (v < 0)
	    *root1 = pow(0.5*(w-v), 1.0/3.0) - (uo3) * pow(2.0 / (w-v), 1.0/3.0) - p / 3.0;
	else
	    *root1 = uo3 * pow(2.0 / (w+v), 1.0/3.0) - pow(0.5*(w+v), 1.0/3.0) - p / 3.0;


	return 1;
    }
  
    if (uo3 >= 0) {//Multiple root detected	  
	*root1 = *root2 = *root3 = pow(v, 1.0 / 3.0) - p / 3.0;
	return 3;
    }

    double muo3 = - uo3;
    double s;
    if (muo3 > 0)
    {
	s = sqrt(muo3);
	if (p > 0) s = -s;
    }
    else
	s = 0;
      
    double scube = s * muo3;
    if (scube == 0)
    {
	*root1 = - p / 3.0;
	return 1;
    }
      
    double t = - v / (scube + scube);
    double k = acos(t) / 3.0;
    double cosk = cos(k);
    *root1 = (s + s) * cosk - p / 3.0;
      
    double sinsqk = 1.0 - cosk * cosk;
    if (sinsqk < 0) return 1;

    double rt3sink = sqrt(3) * sqrt(sinsqk);
    *root2 = s * (-cosk + rt3sink) - p / 3.0;
    *root3 = s * (-cosk - rt3sink) - p / 3.0;

 
    return 3;
}




/**********************************************/
int quadSolve(float C, float B, float A, float* root1, float* root2) {
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
float quarticError(float a, float b, float c, float d,
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
		    errors[root] = pow(abs(value / thirdDeriv), 1.0/3.0);
		else
		    errors[root] = sqrt(sqrt(abs(value)/24));
	    }
	}
	if (max_root < errors[root]) max_root = errors[root];
    }

    return max_root;
}

/**********************************************/
int ferrariQuarticSolve(float a, float b, float c, float d,
			float * root1, float * root2, float * root3, float * root4) {
    float rts[4];
    float worst3[3];
    float qrts[4][3]; /* quartic roots for each cubic root */

    if (d == 0.0) {
	*root1 = 0.0;
	return cubicSolve(a,b,c,root2,root3,root4) + 1;
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
    int n3 = cubicSolve (p, q, r, &v3[0],&v3[1],&v3[2]);
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

	    int n1 = quadSolve(hh, gg, 1.0, &v1[0], &v1[1]);
	    int n2 = quadSolve(h,   g, 1.0, &v2[0], &v2[1]);
	    n4[j3] = n1*2+n2*2;
	    qrts[0][j3] = v1[0];
	    qrts[1][j3] = v1[1];
	    qrts[n1*2+0][j3] = v2[0];
	    qrts[n1*2+1][j3] = v2[1];
	}
	for (j = 0; j < n4[j3]; ++j)
	    rts[j] = qrts[j][j3];

	worst3[j3] = quarticError(a, b, c, d, rts, n4[j3]);
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


/**********************************************************/
/**********************************************************/
int opt_quat_sine_lapack ( double ** x, int NX, int *set_of_directions_x,
	       double ** y, int NY, int *set_of_directions_y,
	       int set_size, double * q, double * rmsd) {

    
    double * x_sub[set_size], * y_sub[set_size];
    int  ctr;
    int  i, j;
 
    double ATA     [4][4] = {{0.0}};
    double prev_ATA[4][4] = {{0.0}};
    double ATA_sum [4][4] = {{0.0}};
    double a[3] = {0.0}, b[3] = {0.0};
    
    int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		       double result[4][4]);
    int construct_ATA (double ATA[4][4], double a[3], double  b[3]);


    if (!set_size) {
	*rmsd = -1;
	return 1;
    }
    memset ( &(q[0]), 0, 4*sizeof(double) );

    /* find the subset */
    ctr = 0;
    for ( ctr =0; ctr < set_size; ctr++ ) {
	x_sub[ctr] =  x[set_of_directions_x[ctr]];
	y_sub[ctr] =  y[set_of_directions_y[ctr]];
    }

    /* B = ATA_sum matrix to diagonalize in order to get the quaternion */
    for ( ctr =0; ctr < set_size; ctr++ ) {
   	for (i=0; i<3; i++ ) {
	    a[i] = y_sub[ctr][i] + x_sub[ctr][i];
	    b[i] = y_sub[ctr][i] - x_sub[ctr][i];
	}
 	construct_ATA (ATA, a, b);
	add_matrices (prev_ATA, ATA, ATA_sum);
	memcpy (prev_ATA[0], ATA_sum[0], 4*4*sizeof(double));
    }
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA_sum[i][j] /= set_size;
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

	float A, B, C, D, E; /* coefficients id the characteristic polynomial - according to sympy */
	A = 1.0;/* will not be used below */
	B = (-a - f - k - p);
	C = (a*f + a*k + a*p - b*e - c*i - d*m + f*k + f*p - g*j - h*n + k*p - l*o);
	D = (- a*f*k - a*f*p + a*g*j + a*h*n - a*k*p + a*l*o + b*e*k + b*e*p - b*g*i
	     - b*h*m - c*e*j + c*f*i + c*i*p - c*l*m - d*e*n + d*f*m - d*i*o + d*k*m
	     - f*k*p + f*l*o + g*j*p - g*l*n - h*j*o + h*k*n) ;
	E =   a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n - b*e*k*p
	    + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m + c*e*j*p - c*e*l*n
	    - c*f*i*p + c*f*l*m + c*h*i*n - c*h*j*m - d*e*j*o + d*e*k*n + d*f*i*o
	    - d*f*k*m - d*g*i*n + d*g*j*m;

	ferrariQuarticSolve(B, C, D, E, &w[0], &w[1], &w[2], &w[3]);

    }
    double min_w = w[0];
    int    min_i = 0;
    for (i=1; i<4; i++) {
	if (w[i] > min_w) continue;
	min_w = w[i];
	min_i = i;
	
    }
    
    for (i=0; i<4; i++) ATA_sum[i][i] -= min_w;
    if (gauss_elim (ATA_sum,  q)) {
	return 1;
    }
    
    /* rmsd is the lowest eignevalue */
    *rmsd = sqrt(min_w);
    /* and the corresponding eigenvector is what we are looking for */
    if (0) {
	/* w contains the eigenvalues */
	printf (" eigenvals without lapack:\n");
	for (i=0; i<4; i++ ) printf ("%8.3lf ", w[i]);
	printf ("\nrmsd: %8.3lf \n", *rmsd);
	printf ("quat:\n");
	for (i=0; i<4; i++ ) printf ("%8.3lf ", q[i]);
	printf ("\n");
	printf ("**********************************\n");
    }
    
    return 0;
}


/**********************************************************/
/**********************************************************/
int opt_quat_lapack ( double ** x, int NX, int *set_of_directions_x,
	       double ** y, int NY, int *set_of_directions_y,
	       int set_size, double * q, double * rmsd) {

    
    double * x_sub[set_size], * y_sub[set_size];
    int  ctr;
    int  i, j;
 
    double ATA     [4][4] = {{0.0}};
    double prev_ATA[4][4] = {{0.0}};
    double ATA_sum [4][4] = {{0.0}};
    double a[3] = {0.0}, b[3] = {0.0};
    
    int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		       double result[4][4]);
    int construct_ATA (double ATA[4][4], double a[3], double  b[3]);

    /* note how we pass the matrix: pointer to the first element in the block */
    void dsyev_ (char * jobz, char *uplo,  int *n,
		  double *A, int * lda, double * w, double * work, int * lwork, int *info);

    if (!set_size) {
	*rmsd = -1;
	return 1;
    }

    memset ( &(q[0]), 0, 4*sizeof(double) );

    
    /* find the subset */
    ctr = 0;
    for ( ctr =0; ctr < set_size; ctr++ ) {
	x_sub[ctr] =  x[set_of_directions_x[ctr]];
	y_sub[ctr] =  y[set_of_directions_y[ctr]];
    }

    /* check: */
    if (0) {
	printf (" Number of vectors to match: %d. \n", set_size);
	for ( ctr =0; ctr < set_size; ctr++ ) {
	    printf ("\t x%1d   %10.4lf  %10.4lf  %10.4lf   ",
		    ctr, x_sub[ctr][0], x_sub[ctr][1], x_sub[ctr][2]);
	    printf ("\t y%1d   %10.4lf  %10.4lf  %10.4lf \n",
		    ctr, y_sub[ctr][0], y_sub[ctr][1], y_sub[ctr][2]);
	}
    }

     
    /* B = ATA_sum matrix to diagonalize in order to get the quaternion */
    for ( ctr =0; ctr < set_size; ctr++ ) {
   	for (i=0; i<3; i++ ) {
	    a[i] = y_sub[ctr][i] + x_sub[ctr][i];
	    b[i] = y_sub[ctr][i] - x_sub[ctr][i];
	}
 	construct_ATA (ATA, a, b);
	add_matrices (prev_ATA, ATA, ATA_sum);
	memcpy (prev_ATA[0], ATA_sum[0], 4*4*sizeof(double));
    }
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA_sum[i][j] /= set_size;
	}
    }
    /* diagonalize ATA_sum - the eigenvector corresponsing to the
       smallest lambda is the quaternion we are looking for; the
       eigenvalue is the rmsd*/
    /* use the nomenclature from dsyev*/
    char jobz= 'V'; /*Compute eigenvalues and eigenvectors.*/
    char uplo= 'U'; /* Upper triangle of A (the matrix we are diagonalizing) is stored; */
    int  n = 4;     /* order and the leading dimension of A */
    int  lda = 4;
    double ** A;
    int  info;
    int  lwork = 200;
    double w [4];
    double work[200];
    
    if ( !( A=dmatrix(4,4)) ) exit (1);
    memcpy (A[0], ATA_sum[0], 4*4*sizeof(double));

   /* note how we pass the matrix: */
    dsyev_ ( &jobz, &uplo,  &n, A[0], &lda, w, work, &lwork, &info);
    if (  ! info) {
	*rmsd = sqrt (w[0]);
	for (i=0; i<4; i++ ) q[i] = A[0][i];
	if (0) {
	    /* w contains the eigenvalues */
	    printf ("**********************************\n");
	    printf ("egeinvalues: \n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", w[i]);
	    printf ("\nrmsd: %8.3lf \n", *rmsd);
	    printf ("quat:\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", q[i]);
	    printf ("\n");
	}
    } else {
	fprintf (stderr, "Error in dsyev().\n");
	return 1;
    }
    
    free_dmatrix(A);
    
    return 0;
}
