/* readin in 2 pdb, presumably roughly aligned structurally, */
/* and output the corresponding sequence alignment */
/* Ivana, Jul 24, 2010 */
# include "nossa.h"

/* the following functions are defined below main()*/

int smith_waterman (int max_i, int max_j, double **similarity,
		    int *map_i2j, int * map_j2i, double * aln_score);
double  two_point_distance (Atom * atom1, Atom *atom2) ;

/*************************************************************************/
/*************************************************************************/
int main (int argc, char * argv[]) {


    Protein protein1, protein2;
    double ** similarity;
    double d, d0 = 10.0, aln_score;

    char pdbname1 [BUFFLEN] = "\0";
    char pdbname2 [BUFFLEN] = "\0";
    int * residue_map_i2j = NULL, * residue_map_j2i = NULL;
    int  resctr1, resctr2;
    int  output_tfm;


    if ( argc < 2 ) {
	printf ( "Usage: %s <pdbname1> <pdbname2> [-tfm].\n", argv[0] );
	exit (1);
    } 
    sprintf ( pdbname1, "%s", argv[1]);
    sprintf ( pdbname2, "%s", argv[2]);

    output_tfm = 0; /* output transformation matrix */

    if (  argc > 3) {
	if ( ! strcmp ( argv[3], "-tfm") ) {
	    output_tfm = 1;
	} else {
	    fprintf (stderr, "Unknonwn argument: %s.\n", argv[3]);
	    exit (1);
	}
    }

    /***************************************/
    /* read in the two structures          */
    if ( read_pdb ( pdbname1, &protein1, '\0'))  exit (1);
    if ( read_pdb ( pdbname2, &protein2, '\0'))  exit (1);
    //printf ("read in %s, length %d\n", pdbname1, protein1.length);
    //printf ("read in %s, length %d\n", pdbname2, protein2.length);

    /**********************************************/
    /* allocate space for the "similarity" matrix */
    similarity = dmatrix (protein1.length, protein2.length);
    if ( ! similarity) exit (1);
    /*******************************************/
    /* allocate space for the maps between the */
    /* sets of Ca's                            */
    if ( residue_map_i2j) free (residue_map_i2j);
    if ( residue_map_j2i) free (residue_map_j2i);
    if ( ! (residue_map_i2j = emalloc (protein1.length*sizeof(int))) ) exit (1); 
    if ( ! (residue_map_j2i = emalloc (protein2.length*sizeof(int))) ) exit (1); 
    /* evaluate the similarity matrix from Ca distances */
    for (resctr1=0; resctr1 < protein1.length; resctr1++) {
	if ( !protein1.sequence[resctr1].Ca) continue;
	for (resctr2=0; resctr2< protein2.length; resctr2++) {
	    if ( !protein2.sequence[resctr2].Ca) continue;
	    d = two_point_distance (protein1.sequence[resctr1].Ca,
				    protein2.sequence[resctr2].Ca);
	    similarity[resctr1][resctr2] = exp(-d/d0);
	}
    }

 

    /* find the maps using SW                 */
    smith_waterman (protein1.length, protein2.length, similarity,
		    residue_map_i2j, residue_map_j2i, &aln_score);
    double fraction = (protein1.length < protein2.length) ?
	aln_score/protein1.length : aln_score/protein2.length;
    printf ("alignment score: %8.3lf %8.3lf \n", aln_score, fraction);
	
 
    return 0;
}


/*************************************************************************/
/*************************************************************************/
int smith_waterman (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char   **direction;
    double gap_opening   = -0.2;
    double gap_extension = -0.1;
    double endgap        =  0.0;
    double far_far_away  = -100.0;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    double F_max;
    int use_endgap = 1;
    int F_max_i, F_max_j;
    int i,j;

     /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;


    
    /* fill the table */
    F_max = far_far_away;
    F_max_i = 0;
    F_max_j = 0;
   
    for (i=0; i<= max_i; i++) {
	for (j=0; j<=max_j; j++) {

	    if ( !i && !j ) {
		F[0][0] = 0;
		direction[i][j] = 'd';
		continue;
	    }
	    
	    if ( i && j ){
		
		if ( direction[i-1][j] == 'i' ) {
		    /*  gap extension  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_extension;		    
		} else {
		    /*  gap opening  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_opening;
		}
		i_sim =  F[i-1][j] + penalty;
		
		if ( direction[i][j-1] =='j' ) {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_extension;		    
		} else {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_opening;		    
		}
		j_sim = F[i][j-1] +  penalty;
       	
		
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] ;
		
		
	    } else if ( j ) {
		
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i][j-1] =='j' ) {
			penalty = gap_extension;
		    } else {
			penalty = gap_opening;
		    }
		}
		j_sim = F[i][j-1] + penalty;
		
		i_sim = diag_sim = far_far_away;

		
	    } else if ( i ) {
		
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i-1][j] == 'i' ) {
			penalty =  gap_extension;
		    } else {
		        penalty =  gap_opening;
		    }
		}
		i_sim = F[i-1][j] + penalty;
		
		j_sim = diag_sim = far_far_away;
		

	    } 

	    max_sim = diag_sim;
	    direction[i][j] = 'd';
	    if ( i_sim > max_sim ){
		max_sim = i_sim;
		direction[i][j] = 'i';
	    }
	    if ( j_sim > max_sim ) {
		max_sim = j_sim;
		direction[i][j] = 'j';
	    }

	    if ( max_sim < 0.0 ) max_sim = 0.0;
	    
	    F[i][j] = max_sim;
	    if ( F_max < max_sim ) {
		/* TODO: tie break here */
		F_max = max_sim;
		F_max_i = i;
		F_max_j = j;
		
	    }
	    
	}
    }
    
    /*retrace from the maximum element*/
    i= max_i;
    for( i= max_i; i>F_max_i; i--) map_i2j[i-1] = far_far_away;;
    for( j= max_j; j>F_max_j; j--) map_j2i[j-1] = far_far_away;;
    i = F_max_i;
    j = F_max_j;
    *aln_score = F[i][j]; 
    while ( i>0 ||  j >0 ) {
	//printf (" %4d  %4d  %8.3f  \n", i, j, F[i][j]);
	switch ( direction[i][j] ) {
	case 'd':
	    //printf ( " %4d  %4d \n",  i, j);
	    map_i2j [i-1] = j-1;
	    map_j2i [j-1] = i-1;
	    i--;
	    j--; 
	    break;
	case 'i':
	    //printf ( " %4d  %4d \n",  i, -1);
	    map_i2j [i-1] = far_far_away;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = far_far_away;
	    j--; 
	    break; 
	default: 
	    fprintf (stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}


/*************************************************************************/
/*************************************************************************/
double  two_point_distance (Atom * atom1, Atom *atom2) {
 
    double xx, yy, zz;;
    double d = 0;
    
    xx = atom1->x - atom2->x; xx *= xx;
    yy = atom1->y - atom2->y; yy *= yy;
    zz = atom1->z - atom2->z; zz *= zz;
    d = xx+yy+zz;
    
    return  sqrt (d);
}
