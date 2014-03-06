# include <struct.h>
# include <float.h>



/**********************************************************/
/**********************************************************/
int opt_quat ( double ** x, int NX, int *set_of_directions_x,
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
