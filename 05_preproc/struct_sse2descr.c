/* /\* */
/* This source code is part of deconSTRUCT, */
/* protein structure database search and backbone alignment application. */
/* Written by Ivana Mihalek, with contributions from Mile Sikic. */
/* Copyright (C) 2012 Ivana Mihalek. */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program. If not, see<http://www.gnu.org/licenses/>. */

/* Contact: ivana.mihalek@gmail.com. */
/* *\/ */

# include "struct.h"




/*************************************************************************************/
int sse2descriptor (Protein *protein, Descr* descr) {


    int resctr, type, prev_type = 0;
    int element_ctr = 0, no_of_sses = 0;
    int point_ctr = 0, first_res_index = 0, last_res_index = 0;
    int sse_id, prev_sse_id;
    double **point = NULL;
    Residue * residue;
 
    if ( ! (point=dmatrix(protein->length,3))) return 1;

    descr->no_of_helices = 0;
    descr->no_of_strands = 0;
    
    /* count the SSEs, and allocate the Description space */
    element_ctr = 0;
    prev_sse_id = 0;
    for (resctr=0; resctr<protein->length; resctr++) {
	residue = protein->sequence+resctr;

	if ( residue->belongs_to_helix) {
	    sse_id = residue->belongs_to_helix;
	}  else if (residue->belongs_to_strand) {
	    sse_id = residue->belongs_to_strand;
	} else {
	    sse_id = 0;
	}

	if ( prev_sse_id && sse_id != prev_sse_id ) element_ctr++;

	prev_sse_id = sse_id;
	
    }
    no_of_sses = element_ctr;
    if ( ! (descr->element = emalloc(no_of_sses*sizeof(SSElement)) )) return 1;


    /* turn SSEs into sets of points to fit onto */
    point_ctr   = 0;
    element_ctr = 0;
    first_res_index = last_res_index = 0;
    sse_id = 0;
    
    for (resctr=0; resctr<protein->length; resctr++) {
	
	residue = protein->sequence+resctr;

	if ( residue->belongs_to_helix) {
	    sse_id = residue->belongs_to_helix;
	    type = HELIX;
	    
	}  else if (residue->belongs_to_strand) {
	    sse_id = residue->belongs_to_strand;
	    type = STRAND;

	} else {
	    sse_id = 0;
	    type = LOOP;
	}

	if ( sse_id != prev_sse_id ) {
	    
	    if ( prev_sse_id  ) {

		last_res_index = resctr-1;
		if (point_ctr >= MIN_CaS_IN_SSE) {
		    /* we didn't process this one */
		    process_sse (protein, prev_type, point, point_ctr,
				 first_res_index, last_res_index,
				 descr->element+element_ctr);
		    element_ctr ++;
		    if (prev_type == HELIX ) {
			descr->no_of_helices ++;
		    } else if (prev_type== STRAND) {
			descr->no_of_strands ++;
		    }
		}
	    }

	    if ( sse_id ) {
		first_res_index = resctr;
		point_ctr = 0;
	    }
	}

	if ( sse_id ) {
	    if ( residue->Ca ) { /* survive cases when Ca not given */
		point[point_ctr][0] = residue->Ca->x;
		point[point_ctr][1] = residue->Ca->y;
		point[point_ctr][2] = residue->Ca->z;
		point_ctr ++;
	    }
	}
	prev_sse_id = sse_id;
	prev_type   = type;
    }

   
    if (prev_sse_id && point_ctr) {
	last_res_index = resctr-1;
	/* we didn't process this one */
	process_sse (protein, prev_type, point, point_ctr, 
		     first_res_index, last_res_index,
		     descr->element+element_ctr);
	if (type == HELIX ) {
	    descr->no_of_helices ++;
	} else if (type== STRAND) {
	    descr->no_of_strands ++;
	}
	element_ctr ++;
    }

    descr->no_of_elements = descr->no_of_helices + descr->no_of_strands;
    descr->no_of_residues = protein->length;
    
    protein->no_helices  = descr->no_of_helices;
    protein->no_strands  = descr->no_of_strands;
    protein->no_of_elements = descr->no_of_elements;

    free_dmatrix(point);
    
    return 0;

}


/********************************************************************/
/********************************************************************/
int process_sse (Protein * protein, int type, double ** point, int number_of_points,
		 int first_res_index, int last_res_index, SSElement * element){

    int point_ctr;
    int number_of_avgd_points;
    double ** point_avg = NULL;
    Residue * residue;
    
    if ( type != HELIX && type !=STRAND) {
	fprintf (stderr, "Error in %s:%d:  attempting a fit "
		 "on something that is neither strand nor helix.\n",
		 __FILE__, __LINE__);
	exit (1);
    }

    /* the fields that are straightforward to fill: */
    element->type    = type;
    
    residue = protein->sequence+first_res_index;
    strncpy (element->begin_id, residue->pdb_id, (PDB_ATOM_RES_NO_LEN+2)*sizeof(char));
    
    residue = protein->sequence+last_res_index;
    strncpy (element->end_id, residue->pdb_id, (PDB_ATOM_RES_NO_LEN+2)*sizeof(char));

    element->chain  = residue->chain;
    element->begin  = first_res_index;
    element->end    =  last_res_index;
    element->length = last_res_index-first_res_index+1;
    
    /* fit a line through the points -- 
       fill in the center and direction fields of the element structure*/
    if (type== STRAND) {
	fit_line (point, number_of_points, element->cm, element->p);
	
    } else {
	number_of_avgd_points = number_of_points-2;
	if ( ! (point_avg=dmatrix(number_of_avgd_points,3))) return 1;
	for (point_ctr = 0; point_ctr < number_of_avgd_points; point_ctr++) {
	    int i;
	    for (i=0; i<3; i++) {
		point_avg[point_ctr][i] = point[point_ctr][i]
		    + point[point_ctr+1][i]
		    + point[point_ctr+2][i];
		point_avg[point_ctr][i] /= 3;
	    }
	    
	}
	fit_line (point_avg, number_of_avgd_points, element->cm, element->p);   
	free_dmatrix (point_avg);   
    }


    return 0;
}




int fit_line (double **point, int no_points, double center[3], double direction[3]) {

    int i, j, point_ctr;
    double I[3][3] = {{0.0}};
    double translated_point[3];
    
    void dsyev_(char * jobz, char * uplo, int* N,
		double * A, int * leading_dim,
		double * eigenvalues,
		double *workspace, int *workspace_size,
		int * retval);
    if (!no_points) {
	fprintf (stderr, "Error in %s:%d: in fit_line() - no points given\n",
		 __FILE__, __LINE__);
	exit (1);
    }
  
    for (i=0; i<3; i++) center[i] = 0;
    for (point_ctr = 0; point_ctr < no_points; point_ctr++) {
	for (i=0; i<3; i++) center[i] += point[point_ctr][i];
    }
    for (i=0; i<3; i++) center[i] /= no_points;
    
    for (point_ctr = 0; point_ctr < no_points; point_ctr++) {

	for (i=0; i<3; i++) {
	    translated_point[i] = point[point_ctr][i] - center[i];
	}
	/* moments of inertia */
	for (i=0; i<3; i++) {  /* modulo = circular permutation */
	    I[i][i] += translated_point[(i+1)%3]*translated_point[(i+1)%3] +
		translated_point[(i+2)%3]*translated_point[(i+2)%3];
	    for ( j=i+1; j<3; j++) { /* offdiag elements */
		I[i][j] -= translated_point[i]*translated_point[j];
	    }
	}
    }
    for (i=0; i<3; i++) {
	for ( j=i+1; j<3; j++) { 
	    I[j][i] =  I[i][j];
	}
    }


    /*****************************************/
    /* diagonalize I[][], pick the direction
       with the smallest moment of inertia,
       and rotate back to the initial frame  */
    /*****************************************/
    char jobz = 'V'; /* find evalues and evectors */
    char uplo = 'L'; /* matrix is stored as lower (fortran convention) */
    int  N = 3; /* the order of matrix */
    int leading_dim = N;
    int retval = 0;
    int workspace_size = 3*N;
    double A[N*N];
    double eigenvalues[N];
    double workspace[3*N];
    double NC_direction[3];
    
    for (i=0; i<3; i++) {
	for ( j=0; j<3; j++) {
	    A[i*3+j] = I[i][j];
	}
    }

    
    dsyev_ ( &jobz, &uplo, &N, A, &leading_dim, eigenvalues,
	     workspace, &workspace_size, &retval);
    
    if ( retval ) {
	fprintf (stderr, "Error in %s:%d: dsyev() returns %d.\n",
		 __FILE__, __LINE__, retval);
	exit (1);
    }
    for (i=0; i<3; i++)  direction[i] = A[i];

    /* sum of squared dist is egeinvalues[0];
       in case we need it */

    /* force the direction to be N-->C */
    for (i=0; i<3; i++) {
	NC_direction[i] = point[no_points-1][i] - point[0][i];
    }

    double dot;
    unnorm_dot (NC_direction,direction, &dot);
    if ( dot<0 ) {
	for (i=0; i<3; i++) direction[i] *= -1;
    }

    return 0;
}

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


