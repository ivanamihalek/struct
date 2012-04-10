# include "struct.h"
# include "structure2sse.h"


int sse2descriptor (Protein *protein, Descr * descr);
int process_sse (int type, double ** points,
		 int number_of_points, SSElement * element);

int main ( int argc, char * argv[]) {

    int retval;
    Protein protein;
    Descr descr;
   
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

    if ( sse2descriptor (&protein, &descr)) {
	fprintf ( stderr, "Error  fitting lines to sse.\n");
	exit (1);
    }


    return 0;
}



/*************************************************************************************/

int sse2descriptor (Protein *protein, Descr* descr) {


    int resctr, type, prev_type = 0;
    int element_ctr, no_of_sses;
    int number_of_points;
    double **points;
    Residue * residue;

    if ( ! (points=dmatrix(protein->length,3))) return 1;
    
    /* count the SSEs, and allocate the Description space */
    element_ctr = 0;
    for (resctr=0; resctr<protein->length; resctr++) {
	residue = protein->sequence+resctr;
	type = 0;
	if ( residue->belongs_to_helix) {
	    type = HELIX;
	}  else if ( residue->belongs_to_strand) {
	    type = STRAND;
	}

	if ( !(type&prev_type) ) element_ctr++;
 	prev_type = type;
    }
    no_of_sses = element_ctr;
    if ( ! (descr->element = emalloc(no_of_sses*sizeof(SSElement)) )) return 1;

    
    /* turn SSEs into sets of points to fit onto */
    number_of_points = 0;
    for (resctr=0; resctr<protein->length; resctr++) {
	
	residue = protein->sequence+resctr;

	type = 0;
	if ( residue->belongs_to_helix) {
	    type = HELIX;
	}  else if (residue->belongs_to_strand) {
	    type = STRAND;
	}

	if ( ! (type&prev_type) ) {
	    if (prev_type) {
		/* we didn't process this one */
		process_sse (type, points, number_of_points, descr->element+element_ctr);
	    }
	    element_ctr ++;
	    number_of_points = 0;
	}  
        number_of_points ++;
	prev_type = type;
    }

   
    if (prev_type&(HELIX|STRAND) ) {
	/* we didn't process this one */
	process_sse (prev_type, points, number_of_points, descr->element+element_ctr);
    }


    free_dmatrix(points);
    
    return 0;

}

int process_sse (int type, double ** points,
		 int number_of_points, SSElement * element){

    
    /* fit a line through the points */
    

    /* store as Description object */

    return 0;
}



# if 0

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


