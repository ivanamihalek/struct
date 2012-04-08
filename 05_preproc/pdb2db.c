# include "struct.h"
# include "structure2sse.h"


int main ( int argc, char * argv[]) {

    int retval;
    //int seq_length = 0;
    Protein protein;
   
    if ( argc < 3 ) {
	fprintf ( stderr,
		  "Usage: %s <pdb file> [<chain>].\n",
		  argv[0]);
	exit (1);
    }

    retval =  read_pdb(argv[1], argv[2][0], &protein);
    if ( retval ) {
	fprintf ( stderr, "Error reading %s, retval %d.\n",
		  argv[1], retval);
	exit (1);
    }


    if ( structure2sse (&protein)) {
	fprintf ( stderr, "Error  finding SSEs..\n");
	exit (1);
    }

# if 0
    if ( sse2db (&protein)) {
	fprintf ( stderr, "Error  fitting lines to sse.\n");
	exit (1);
    }

 # endif   
    

    
    return 0;
}
/*************************************************************************************/

# if 0
int sse2db (Protein *protein) {

    /* turn SSEs into sets of points to fit onto */

    /* foreach SSE, do the fit */

    /* store as reduced representation */


    return 0;

}


int fit_line (double *points, int no_points, double center[3], double direction[3]) {

    /* find the "moments of inertia" */
    double I[3][3] = {{0.0}};
    
    void dsyev_(char * jobz, char * uplo, int* N,
		double * A, int * leading_dim,
		double * eigenvalues,
		double *workspace, int *workspace_size,
		int * retval);
    
    for (int p = bgn; p <= end; p++) {
	double my_point[3] = {points[p].x - pos.x,
			      points[p].y - pos.y,
			      points[p].z - pos.z};
	for (int x = 0; x < 3; x++) {  /* modulo = circular permutation */
	    I[x][x] += my_point[(x+1)%3]*my_point[(x+1)%3] +
		my_point[(x+2)%3]*my_point[(x+2)%3];
	    for (int y = x+1; y < 3; y++) { /* offdiag elements */
		I[x][y] -= my_point[x]*my_point[y];
	    }
	}
    }
    for (int x = 0; x < 3; x++) { 
	for (int y = x+1; y < 3; y++) {
	    I[y][x] =  I[x][y];
	}
    }

    /*****************************************/
    /* diagonalize I[][], pick the direction
       with the smallest moment of inertia,
       and rotate back to the initial frame  */
    /*****************************************/
    char jobz = 'V'; /* find evalues and evectors */
    char uplo = 'L'; /* amtrix is stored as lower (fortran convention) */
    int  N = 3; /* the order of matrix */
    int leading_dim = N;
    int retval;
    double A[N*N];
    double eigenvalues[N];
    double workspace[3*N];
    int workspace_size = 3*N;

    for (int x = 0; x < 3; x++) {
	for (int y=0; y < 3; y++) {
	    A[x*3+y] = I[x][y];
	}
    }
   
    dsyev_ ( &jobz, &uplo, &N, A,  &leading_dim, eigenvalues,
	     workspace, &workspace_size, &retval);
    if ( retval ) {
	cerr << "Dsyev  error: " << retval << "\n";
	throw(MLException("eigen decompfailed"));
    }


}
# endif

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
# if 0
/* junkyard */
    int  resctr;
    int ord_number;
    for (resctr=0; resctr<protein.length; resctr++) {
        ord_number = -1;
	char sse_code = 'C';
	if ( protein.sequence[resctr].belongs_to_helix) {
	    sse_code = 'H';
            ord_number = protein.sequence[resctr].belongs_to_helix;
	} else if ( protein.sequence[resctr].belongs_to_strand) {
	    sse_code = 'E';
            ord_number = protein.sequence[resctr].belongs_to_strand;
	}
        	
	printf ("  %s  %c  %d\n", protein.sequence[resctr].pdb_id, sse_code, ord_number);

	
    }
    


# endif


