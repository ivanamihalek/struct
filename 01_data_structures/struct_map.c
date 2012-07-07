/*
This source code is part of deconSTRUCT,
protein structure database search and backbone alignment application.
Written by Ivana Mihalek, with contributions from Mile Sikic.
Copyright (C) 2012 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# include "struct.h"


/********************************/
int initialize_map (Map *map, int NX, int NY ) {

    memset (map, 0, sizeof(Map));
    map->x2y = emalloc (NX*sizeof(int));
    map->y2x = emalloc (NY*sizeof(int));
    map->x2y_size = NX;
    map->y2x_size = NY;
    map->cosine  = dmatrix (NX, NY);
    map->sse_pair_score   = dmatrix (NX, NY);
    if ( !map->x2y ||  !map->y2x ||
	 !map->cosine || !map->sse_pair_score ) {
	return 1;
    }
    map->submatch_best = NULL;


    return 0;
}
/********************************/
int clear_map (Map *map) {

    int NX = map->x2y_size;
    int NY = map->y2x_size;

    int  *hide1[2];;
    double**hide2[2];

    hide1[0] = map->x2y;    hide1[1] = map->y2x;
    hide2[0] = map->cosine; hide2[1] = map->sse_pair_score;

    memset (map, 0, sizeof(Map));
    
    map->x2y            = hide1[0];
    map->y2x            = hide1[1];
    map->cosine         = hide2[0];
    map->sse_pair_score = hide2[1];
    
    
    memset (map->x2y, 0, NX*sizeof(int));
    memset (map->y2x, 0, NY*sizeof(int));
    memset (map->cosine[0], 0,  NX*NY*sizeof(double));
    memset (map->sse_pair_score[0], 0,  NX*NY*sizeof(double));

    map->x2y_size = NX;
    map->y2x_size = NY;

    return 0;
}



/********************************/
int copy_map (Map *to, Map* from) {

    if ( from->x2y_size !=  to->x2y_size
	 ||  from->y2x_size !=  to->y2x_size) {
	fprintf ( stderr, "%s:%d:  two maps to be copied "
		  "must be of the same size.\n", __FILE__, __LINE__);
	exit (1);
    }
    
    int NX = from->x2y_size;
    int NY = from->y2x_size;

    to->size    = from->size;
    to->matches = from->matches;
    
    memcpy (to->q, from->q, 4*sizeof(double));
    memcpy (to->T, from->T, 3*sizeof(double));
    
    memcpy (to->x2y, from->x2y, NX*sizeof(int));
    memcpy (to->y2x, from->y2x, NY*sizeof(int));

    to->avg_length_mismatch = from->avg_length_mismatch;
    to->rmsd = from->rmsd;

    
    memcpy (&(to->cosine[0]), &(from->cosine[0]),  NX*NY*sizeof(double));    
    memcpy (&(to->sse_pair_score[0]), &(from->sse_pair_score[0]),  NX*NY*sizeof(double));

  
    to->F = from->F;
    to->avg = from->avg;
    to->avg_sq = from->avg_sq;

    to->z_score = from->z_score;
    to->assigned_score = from->assigned_score;

    if ( from->x2y_residue_level ) {

	fprintf (stderr, "copying res level info in a map not implemented\n");
	exit (1);
    }
   

    if ( from->submatch_best ) {

	fprintf (stderr, "copying submatch info in a map not implemented\n");
	exit (1);
    }

    if ( from->filename[0] ) {
	sprintf (to->filename, "%s", from->filename);
    }

    
    return 0;
}






/********************************/
int free_map (Map *map) {

    if ( !map) return 0;
    
    free (map->x2y);
    free (map->y2x);
    map->x2y_size = 0;
    map->y2x_size = 0;
    
    if (map->x2y_residue_level)
	free (map->x2y_residue_level);
    if (map->y2x_residue_level)
	free (map->y2x_residue_level);
    map->x2y_residue_l_size = 0;
    map->y2x_residue_l_size = 0;
 
    
    free_dmatrix (map->cosine);
    free_dmatrix (map->sse_pair_score); 
    if (map->submatch_best) {
	free (map->submatch_best );
    }


    return 0;
}
/*********************************************************************/  
int list_alloc (List_of_maps * list, int NX, int NY, int fake) {
    
    int map_ctr;

    list->no_maps_allocated      = MAP_MAX*9;
    list->no_maps_used           = 0;
    list->best_array_allocated =   list->no_maps_allocated;
    list->best_array_used      = 0;
    
    list->map = emalloc (list->no_maps_allocated*sizeof(Map) );/* TODO is this enough? */
    if ( !list->map) return 1;
    
    list->map_best    = emalloc (list->best_array_allocated*sizeof(int));
    if (!list->map_best) return 1;
    
    list->NX_allocated = NX;
    list->NY_allocated = NY;

    if (fake) { /* fake init - these will actually be used as pointers
		   to containers in other maps */
    } else {
	
	
	for ( map_ctr= 0; map_ctr<list->no_maps_allocated; map_ctr++) {
	    if ( initialize_map(list->map+map_ctr, NX, NY) ) return 1;
	}
    }

    return 0;
}


int list_report (List_of_maps * list) {

    int map_ctr, best_ctr;

    printf ("==========================================\n");
    printf ("list 0x%x\n", list);
    printf ("maps allocated:   %4d   maps used : %4d   best used : %4d  \n",
	    list->no_maps_allocated, list->no_maps_used, list->best_array_used);


    if ( list->best_array_used ) {
	
	for ( best_ctr= 0; best_ctr<list->no_maps_used; best_ctr++) {

	    map_ctr = list->map_best[best_ctr];
	    
	    printf ("\n  %3d: %3d   %8.3lf %8.3lf %8.3lf %8.3lf \n",
		    best_ctr, map_ctr, 
		    list->map[map_ctr].q[0],
		    list->map[map_ctr].q[1],
		    list->map[map_ctr].q[2],
		    list->map[map_ctr].q[3]);

	    print_map (stdout, list->map+map_ctr, NULL, NULL, NULL, NULL, 2);
	}
 

    } else {
    
	for ( map_ctr= 0; map_ctr<list->no_maps_used; map_ctr++) {
	
	    printf ("  %3d   %8.3lf %8.3lf %8.3lf %8.3lf \n", map_ctr,
		    list->map[map_ctr].q[0],
		    list->map[map_ctr].q[1],
		    list->map[map_ctr].q[2],
		    list->map[map_ctr].q[3]);

	    print_map (stdout, list->map+map_ctr, NULL, NULL, NULL, NULL, 2);
	}
    }

 
    
    return 0;

}

    
/*********************************************************************/  
int list_shutdown (List_of_maps * list, int fake) {
    
    if ( !list || !list->map ) return 0;
    int map_ctr;

    if ( fake ) {
	for ( map_ctr= 0; map_ctr< list->no_maps_used; map_ctr++) {
	    Map *current_map = list->map + map_ctr;
	    if (current_map->x2y_residue_level ) free (current_map->x2y_residue_level);
	    if (current_map->y2x_residue_level ) free (current_map->y2x_residue_level);
	}
    } else {
	for ( map_ctr= 0; map_ctr< list->no_maps_allocated; map_ctr++) {
	    if ( free_map(list->map+map_ctr) ) return 1;
	}
    }
    free (list->map_best);
    free (list->map);
    list->map = NULL;
	
    return 0;
}

/********************************/
/********************************/
/********************************/

int store_sse_pair_score ( Representation *X_rep,  Representation *Y_rep, 
		  double **R, double alpha,  Map *map ){
    
    int i,j, i_is_helix;
    int NX=X_rep->N_full, NY=Y_rep->N_full;
    double cosine;
    double alpha_sq = alpha*alpha;
    double ** x_full_rotated = NULL, **y_full = Y_rep->full;
    double length_dif, exponent;

    if ( ! (x_full_rotated=dmatrix (NX,3) ) ) return 1;
    
    rotate (x_full_rotated, X_rep->N_full, R, X_rep->full);
    
    for (i=0; i < NX; i++ ) {
	i_is_helix =  ( X_rep->full_type[i] == HELIX);
	for (j=0; j < NY; j++) {

	    if ( X_rep->full_type[i] == Y_rep->full_type[j] ) {
		unnorm_dot (x_full_rotated[i], y_full[j], &cosine);
		exponent = 2*(1.0-cosine)/alpha_sq;
		if ( exponent < 10 ) {
		    map->sse_pair_score[i][j]  = exp (-exponent);
		    /* are we using the info about the length? */
		    if (options.use_length && cosine > 0.0) {
			length_dif = abs (X_rep->length[i] - Y_rep->length[j]);
			if ( i_is_helix ) {
			    map->sse_pair_score[i][j]  *=  1/(1.0+length_dif/options.H_length_mismatch_tol);
			} else {
			    map->sse_pair_score[i][j]  *=  1/(1.0+length_dif/options.S_length_mismatch_tol);
			}
		    }
		} else {
		    map->sse_pair_score[i][j] = 0.0;
		}
		map->cosine[i][j] = cosine;
	    } else {
		map->sse_pair_score[i][j]  = options.far_far_away;
		map->cosine[i][j] = -0.9999;
	    }
	}
    }
    free_dmatrix (x_full_rotated);
    
    return 0;
}


/***************************************/
/***************************************/
/***************************************/
int construct_translation_vecs ( Representation *X_rep,  Representation *Y_rep,
				 Map *map ){

    int i, j, k, map_size;
    int NX=X_rep->N_full, NY=Y_rep->N_full;
    double norm_x, norm_y, aux;
 
    /********************************************/
    /* find origin (geom center of mapped SSEs) */
    for (k=0; k<3; k++) {
	X_rep->origin[k] = 0;
	Y_rep->origin[k] = 0;
    }

    map_size = 0;
    for (i=0; i<NX; i++ ) {
	
	j =  map->x2y[i];
	if ( j < 0 ) continue;

	for (k=0; k<3; k++) {
	    X_rep->origin[k] += X_rep->cm[i][k];
	    Y_rep->origin[k] += Y_rep->cm[j][k];
	}
	map_size++;
    }

    for (k=0; k<3; k++) {
	X_rep->origin[k] /= map_size;
	Y_rep->origin[k] /= map_size;
    }

    
    /*******************************************/
    /* find translation vectors                */
    memset (X_rep->translation[0], 0, 3*NX*sizeof(double));
    memset (Y_rep->translation[0], 0, 3*NY*sizeof(double));
    
    for (i=0; i<NX; i++ ) {
	j =  map->x2y[i];
	if ( j < 0 ) continue;
	
	norm_x = norm_y = 0.0;
	for (k=0; k<3; k++) {
	    X_rep->translation[i][k] = aux = X_rep->cm[i][k] -  X_rep->origin[k];
	    norm_x += aux*aux;
	    Y_rep->translation[j][k] = aux = Y_rep->cm[j][k] -  Y_rep->origin[k];
	    norm_y += aux*aux;
	}
	norm_x = sqrt (norm_x);
	norm_y = sqrt (norm_y);

	X_rep->transl_norm[i] = norm_x;
	Y_rep->transl_norm[j] = norm_y;
    }
    
    return 0;
    
}
/***************************************/
/***************************************/
/***************************************/
int find_map ( Penalty_parametrization * penalty_params,
	       Representation *X_rep,  Representation *Y_rep,
	       double ** R,
	       double alpha, double * F_effective,  Map *map,
	       int *anchor_x, int *anchor_y, int anchor_size) {
    
    int i, j, map_size;
    int NX=X_rep->N_full, NY=Y_rep->N_full;
    double aln_score;
    
    /* "sse_pair_score" is the matrix of values (exp weight)*(-1 for SSE mismatch) */
    store_sse_pair_score (X_rep, Y_rep, R,  alpha, map);

    /* if the anchors are given, somebody insists
       they should be considered as already aligned */
    if (penalty_params->custom_gap_penalty_x)
	memset (penalty_params->custom_gap_penalty_x, 0, NX*sizeof(double) );
    if (penalty_params->custom_gap_penalty_y)
	memset (penalty_params->custom_gap_penalty_y, 0, NY*sizeof(double) );
    
    if (anchor_x && anchor_y) {
	 int a;
         for (a=0; a<anchor_size; a++){
	     
	      i = anchor_x[a];
	      if (penalty_params->custom_gap_penalty_x)
		  penalty_params->custom_gap_penalty_x[i] = options.far_far_away;
	      for (j=0; j<NY; j++) {
		   if ( j ==  anchor_y[a] )  continue;
		   map->sse_pair_score[i][j]  =  options.far_far_away;
		   map->cosine[i][j] = -0.9999;
	      }
	      
	      j = anchor_y[a];
	      if (penalty_params->custom_gap_penalty_y)
		  penalty_params->custom_gap_penalty_y[j] = options.far_far_away;
	      for (i=0; i<NX; i++) {
		   if ( i ==  anchor_x[a] )  continue;
		   map->sse_pair_score [i][j] =  options.far_far_away;
		   map->cosine[i][j] = -0.9999;
	      }

	 }
    }
    
    /* dynamic programming using the sse_pair_score */
    if (options.current_algorithm == SEQUENTIAL) {
        smith_waterman (penalty_params, NX, NY, map->sse_pair_score, map->x2y, map->y2x, &aln_score);
    } else if  (options.current_algorithm == OUT_OF_ORDER) {
         hungarian_alignment (NX, NY, map->sse_pair_score, map->x2y, map->y2x, &aln_score);
    } else {
        printf("Wrong algorithm type\n");
        return 1;
    }
   
    map_assigned_score (X_rep, map);

    map -> avg_length_mismatch = 0;
    map_size = 0;
    for (i=0; i<NX; i++) {
	 j = map->x2y[i];
	 if ( j<0) continue;
	 
	 if ( map->cosine[i][j] < options.far_away_cosine) {
	     map->x2y[i] =  options.far_far_away; 
	     map->y2x[j] =  options.far_far_away; /*everything crashes; investigate later */
	 } else {
	     map_size++;
	     map->avg_length_mismatch += fabs(X_rep->length[i] - Y_rep->length[j]);
	     
	 }    
    }

    map->avg_length_mismatch /= map_size;
    *F_effective = -map_size;
    
    return 0;
}

/************************************/
/************************************/
int map_assigned_score ( Representation *X_rep,  Map* map) {

    int i,j;
    int NX = X_rep->N_full;
    
    map->size = 0;
    map->assigned_score = 0;
    for (i=0; i < NX; i++ ) {
	j =  map->x2y[i];
	if ( j <  0 )  continue;
	
	map->assigned_score += map->sse_pair_score[i][j];
    }
    
    /* what's the difference btw this and map->size?*/
    map->matches = match_length (NX, map->x2y);
    return 0;    
}

/************************************/
/************************************/
double find_map_overlap (Map *map1, Map *map2){

    int i, j1, j2;
    double **a= map1->sse_pair_score, **b= map2->sse_pair_score;
    double overlap, norm1, norm2;
    
    if ( map1->x2y_size != map2->x2y_size ) {
	fprintf (stderr, "%s:%d: map size mismatch %d, %d (?).\n",
		 __FILE__, __LINE__, map1->x2y_size, map2->x2y_size);
	exit (1);
    }

    overlap = 0;
    norm1   = 0;
    norm2   = 0;
    for (i=0; i<map1->x2y_size; i++) {
	j1 = map1->x2y[i];
	j2 = map2->x2y[i];
	if ( j1 >=0 && j2 >= 0 && j1==j2) {
	    overlap += a[i][j1]*b[i][j2];
	}
	if ( j1 >= 0 ) {
	    norm1 += a[i][j1]*a[i][j1];
	}
	if ( j2 >= 0 ) {
	    norm2 += b[i][j2]*b[i][j2];
	}
    }

    if ( !norm1 || !norm2 ) {

	fprintf (stderr, "%s:%d: empty map (?).\n", __FILE__, __LINE__);
	exit (1);
    }

    overlap /= sqrt(norm1*norm2);
    
    
    return overlap; 
}




/************************************/
/************************************/
int find_uniq_maps (List_of_maps  *list1, List_of_maps *list2, List_of_maps  *list_uniq) {


    int best_map_ctr, uniq_map_ctr, list_ctr;
    int overlap_found;
    int start[2] = {0,0};
    Map *current_map;
    List_of_maps  *list[2] = {list1, list2};

    uniq_map_ctr      = 0;

    if ( list1) {
	memcpy (list_uniq->map+0, list1->map + list1->map_best[0], sizeof(Map));
	start[0] = 1;
	start[1] = 0;
	
    } else if (list2) {
	memcpy (list_uniq->map+0, list2->map + list2->map_best[0], sizeof(Map));
	start[0] = 0;
	start[1] = 1;
	
    } else {
	fprintf (stderr, "Error in %s:%d: Both lists empty in find_uniq_maps.\n",
		 __FILE__, __LINE__);
	exit (1);
    }
    
    list_uniq->no_maps_used = 1;

    for (list_ctr=0; list_ctr<2; list_ctr++) {
	
	if ( ! list[list_ctr] ) continue;

	
	for ( best_map_ctr= start[list_ctr];  best_map_ctr<list[list_ctr]->no_maps_used;  best_map_ctr++) {

	    if (list[list_ctr]->map_best[best_map_ctr] < 0) break;
	    if (list[list_ctr]->map_best[best_map_ctr] >= list[list_ctr]->best_array_used) break;
	    current_map = list[list_ctr]->map + list[list_ctr]->map_best[best_map_ctr];

	    overlap_found = 0;
	    /* for all unique maps: does current_map overlap? */
	    for ( uniq_map_ctr= 0;  uniq_map_ctr<list_uniq->no_maps_used;  uniq_map_ctr++) {
		/* if current overlaps uniq, store that info  */
		double map_overlap = find_map_overlap (list_uniq->map + uniq_map_ctr, current_map);

		
		if ( map_overlap >= MAP_SIM_THRESHOLD) {
		    overlap_found = 1;
		    break;
		} 
	    }

	    if (!overlap_found) {
		memcpy (list_uniq->map+list_uniq->no_maps_used, current_map, sizeof(Map));
		list_uniq->no_maps_used ++;
	    }
	}

       
    }

    
    for ( uniq_map_ctr= 0;  uniq_map_ctr<list_uniq->no_maps_used;  uniq_map_ctr++) {
	/* these should already be sorted */
	list_uniq->map_best[uniq_map_ctr] = uniq_map_ctr;
    }
    list_uniq->best_array_used = uniq_map_ctr;
    list_uniq->map_best[uniq_map_ctr] = -1;
    
    return 0;
}


/************************************/
/************************************/

/**
 * Function that calculate alignment, and maps structures using Hungarian algorithm
 * @param NX
 * @param NY
 * @param similarity
 * @param x2y
 * @param y2x
 * @param alignment
 * @return 
 */

int hungarian_alignment (int NX, int NY, double **similarity, int * x2y, int * y2x, double * alignment ) {
    int i,j;
    *alignment  = 0;
    for (i =0; i < NX; ++i) x2y[i] = options.far_far_away;
    for (i =0; i < NY; ++i) y2x[i] = options.far_far_away;
     
    
    int multiplier = 1000; // precision level for conversion of double to int
    int **scoring_matrix;
    scoring_matrix = intmatrix(NX, NY);
    
    similarity_to_scoring(NX, NY, multiplier, similarity, scoring_matrix );
    
    
    hungarian_problem_t p;
    
    int matrix_size = hungarian_init(&p, scoring_matrix, NX,NY, HUNGARIAN_MODE_MAXIMIZE_UTIL);
    if (matrix_size < 0) printf("wrong matrix size for Hungarian algorithm\n");

    hungarian_solve(&p);

    // checking if the number of rows is greater than number of columns. Note: Hungarian 
    
    
    // if (NX >= NY) {
    for (i = 0; i < NX; ++i) {
        for (j = 0; j < NY; ++j) {
            if (p.assignment[i][j] && similarity[i][j] > 0) {
                x2y[i] = j;
                y2x[j] = i;
                *alignment += similarity[i][j];
            }
        }
    }
    //} 
        
  
    /* free used memory */
    hungarian_free(&p);
    free_imatrix(scoring_matrix);
    
    return 0;
    
}

/**
 * Function that converts similarity matrix of doubles to scoring matrix of integers
 * Scoring matrix is necessary to the Hungarian algorithm 
 * @param m input double matrix
 * @param rows
 * @param cols
 * @return matrix of integers
 */

void similarity_to_scoring(int NX, int NY, int multiplier, double** similarity, int ** hungarian_alignment ) {
  int i,j;
  
  for(i=0;i<NX;i++) {
      for(j=0;j<NY;j++) hungarian_alignment[i][j] = similarity[i][j] * multiplier;
  }
}



/************************************/
/************************************/
int smith_waterman (Penalty_parametrization *params, int max_i, int max_j, double **similarity,
		    int *map_i2j, int * map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char ** direction;
    double *custom_gap_penalty_x =  params->custom_gap_penalty_x;
    double *custom_gap_penalty_y =  params->custom_gap_penalty_y;
    double gap_opening   = params->gap_opening;
    double gap_extension = params->gap_extension;
    double endgap        = params->endgap;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    double F_max;

    int use_endgap = params->endgap_special_treatment;
    int F_max_i, F_max_j;
    int i,j;

     /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;

    for (i =0; i < max_i; ++i) map_i2j[i] = options.far_far_away;
    for (j =0; j < max_j; ++j) map_j2i[j] = options.far_far_away;
    
    
    /* fill the table */
    F_max = options.far_far_away;
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
		    if ( custom_gap_penalty_x && custom_gap_penalty_x[i-1] < 0 ) {
			penalty = custom_gap_penalty_x[i-1];
		    } else {
			penalty = (use_endgap&&j==max_j) ? endgap : gap_extension;
		    }
		} else {
		    /*  gap opening  */
		    if ( custom_gap_penalty_x && custom_gap_penalty_x[i-1] < 0 ) {
			penalty = custom_gap_penalty_x[i-1];
		    } else {
			penalty = (use_endgap&&j==max_j) ? endgap : gap_opening;
		    }
		}
		i_sim =  F[i-1][j] + penalty;
		
		if ( direction[i][j-1] == 'j' ) {
		    if ( custom_gap_penalty_y && custom_gap_penalty_y[j-1] < 0 ) {
			penalty = custom_gap_penalty_y[j-1];
		    } else {
			penalty = (use_endgap&&i==max_i) ? endgap : gap_extension;
		    }
		} else {
		    if ( custom_gap_penalty_y && custom_gap_penalty_y[j-1] < 0 ) {
			penalty = custom_gap_penalty_y[j-1];
		    } else {
			penalty = (use_endgap&&i==max_i) ? endgap : gap_opening;
		    }
		}
		j_sim = F[i][j-1] +  penalty;
       	
		
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] ;
		

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
		
	    } else if ( j ) {
		
		if ( custom_gap_penalty_y && custom_gap_penalty_y[j-1] < 0 ) {
		    penalty = custom_gap_penalty_y[j-1];
		} else {
		    if ( use_endgap) {
			penalty = endgap;
		    } else {
			if ( direction[i][j-1] =='j' ) {
			    penalty = gap_extension;
			} else {
			    penalty = gap_opening;
			}
		    }
		}
		j_sim = F[i][j-1] + penalty;
		max_sim = j_sim;
		direction[i][j] = 'j';

		
	    } else if ( i ) {
		
		if ( custom_gap_penalty_x && custom_gap_penalty_x[i-1] < 0 ) {
		    penalty = custom_gap_penalty_x[i-1];
		} else {
		    if ( use_endgap) {
			penalty = endgap;
		    } else {
			if ( direction[i-1][j] == 'i' ) {
			    penalty =  gap_extension;
			} else {
			    penalty =  gap_opening;
			}
		    }
		}
		i_sim = F[i-1][j] + penalty;
		max_sim = i_sim;
		direction[i][j] = 'i';
		

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
    for( i= max_i; i>F_max_i; i--) map_i2j[i-1] = options.far_far_away;;
    for( j= max_j; j>F_max_j; j--) map_j2i[j-1] = options.far_far_away;;
    i = F_max_i;
    j = F_max_j;
    *aln_score = F[i][j]; 
    while ( i>0 ||  j >0 ) {
	//printf (" %4d  %4d  %8.3f  \n", i, j, F[i][j]);
	if ( i<0 || j<0 ) {
	    fprintf ( stderr, "Retracing error.\n");
	}
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
	    map_i2j [i-1] = options.far_far_away;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = options.far_far_away;
	    j--; 
	    break; 
	default: 
	    fprintf ( stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}
/************************************/
/************************************/
int needleman_wunsch (int max_i, int max_j, double **similarity,
		      int *map_i2j, int *map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char ** direction;
    double gap_opening   = options.gap_open;
    double gap_extension = options.gap_extend;
    double endgap        = options.endgap;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    int use_endgap = options.use_endgap;
    int i,j;

    /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;


    /* fill the table */
    for (i=0; i<= max_i; i++) {
	for (j=0; j<=max_j; j++) {

	    if ( !i && !j ) { /* upper left corner */
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
		
		i_sim = diag_sim = options.far_far_away;

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
		
		j_sim = diag_sim = options.far_far_away;
		
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

	    F[i][j] = max_sim;
	    
	}
    }
    *aln_score = F[max_i][max_j];
    
    /*retrace*/
    i = max_i;
    j = max_j;
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
	    map_i2j [i-1] = options.far_far_away;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = options.far_far_away;
	    j--; 
	    break; 
	default: 
	    fprintf ( stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}

/************************************************************/
/************************************************************/

int match_length (int N, int * x2y) {
    int i, length = 0;
    for (i=0; i<N; i++ ) {
	if (x2y[i] >= 0 ) length ++;
    }
    return length;
}
/********************************/
/********************************/
/********************************/
int map_complementarity (Map *map1, Map *map2,  double *z) {
    
    int i, j, size = 0;
    int * map[2];
    double **array;
    double **sse_pair_score[2];
    double complementarity = 0.0, z_c,avg, avg_sq, stdev;
    double upper, lower;
    int ctr, length;

    size =  map1->x2y_size +  map2->y2x_size;
    array = dmatrix(2, size);
    if (! array ) return 1;
    
    upper = -100;
    lower = 100;
    ctr = 0;

    map[0]  =  map1->x2y;
    map[1]  =  map2->x2y;
    sse_pair_score[0] = map1->sse_pair_score;
    sse_pair_score[1] = map2->sse_pair_score;
    
	
    length = 0;
    for (i=0; i < map1->x2y_size ; i++ ) {
	if ( map[0][i] < 0 &&  map[1][i] < 0 ) continue;
	
	for (ctr=0; ctr<2; ctr++ ) {
	    if ( (j=map[ctr][i]) >= 0 ) {
		array[ctr][length] = sse_pair_score[ctr][i][j];
		if (upper < array[ctr][length] ) {
		    upper =  array[ctr][length];
		} else if ( lower > array[ctr][length] ) {
		    lower =  array[ctr][length];
		}
	    } else {
		array[ctr][length] = -10; 
	    }
	}
	length++;	
    }
    for (i=0; i<length; i++ ) {
	for (ctr=0; ctr<2; ctr++ ) {
	    if ( array[ctr][i] < -1 ) {
		array[ctr][i] = lower;
	    }
	}
    }

    /* printf (">>>>>>  l, u %6.2lf %6.2lf\n", lower, upper); */
   
    if ( length > 3) {
	double diff;
	double avg_1, avg_2, avg_3, avg_4;
	double u = upper, l = lower;
	double interval = u-l;
	int N = length;
	
	complementarity = 0.0;
	for (i=0; i < length; i++ ) {
	    /* printf (">>>   %2d    %6.2lf %6.2lf\n", i, */
	    /* 	    array[0][i], array[1][i]); */
	    diff = array[0][i] - array[1][i] ;
	    complementarity += diff*diff;
	}
	
	avg_1 = (u*u - l*l)/2/interval;
	avg_2 = (u*u*u - l*l*l)/3/interval;
	avg_3 = (u*u*u*u - l*l*l*l)/4/interval;
	avg_4 = (u*u*u*u*u - l*l*l*l*l)/5/interval;
   
	avg = 2*N*( avg_2 - avg_1*avg_1);
	avg_sq = N*(N-1)*( 4*avg_2*avg_2 - 8*avg_2*avg_1*avg_1 +
			   4*avg_1*avg_1*avg_1*avg_1)
	    + N*(2*avg_4 - 8*avg_3*avg_1 + 6*avg_2*avg_2);
	stdev = sqrt (avg_sq - avg*avg);
	z_c = (complementarity-avg)/stdev;
	/* printf (" complementarity =   %6.2lf   z =   %6.2lf \n", */
        	  /*   complementarity/length, z_c); */
	*z = z_c;
    } else {
	*z = -10;
    }

    
    //if ( length > 3) kendl1(array1, array2, length, tau, z, prob);

    free_dmatrix(array);
    
    return  0;
    
    
}
/********************************/
/********************************/
/********************************/
int map_consistence (  int NX, int NY, Map * combined_map, Map *map1, Map *map2,
		       double *total_ptr,  double *gap_score) {
    int i,j;
    double val1, val2;
    double total = 0;
    double aln_score;

    if (!NX) NX = map1->x2y_size; /* TODO: rename */
    if (!NY) NY = map1->y2x_size; /* TODO: rename */
    
    for (i=0; i<NX; i++) {
	for (j=0; j<NY; j++) {
	    val1 =  map1->sse_pair_score[i][j];
	    val2 =  map2->sse_pair_score[i][j];

	    if ( val1 > val2) {
		combined_map->sse_pair_score[i][j]  = val1;
		combined_map->cosine[i][j] = map1->cosine[i][j];
	    } else {
		combined_map->sse_pair_score[i][j]  = val2;
		combined_map->cosine[i][j] = map2->cosine[i][j];
	    }
	}
    }
      

    /* Needleman on combined sse_pair_score */
    needleman_wunsch (NX, NY, combined_map->sse_pair_score,
		      combined_map->x2y,  combined_map->y2x, &aln_score );
    /* how do lengths compare? */
    combined_map->matches = match_length (NX, combined_map->x2y);

    if ( ! combined_map->matches ) {
	combined_map->assigned_score = 0.0;
    } else {
	for (i=0; i<NX; i++) {
	    j = combined_map->x2y[i];
	    if (j<0) continue;
	    total += combined_map->sse_pair_score[i][j];
	}
	combined_map->assigned_score = total;
    }
    /* pass NULL for total_ptr  if I don't need it */
    if (total_ptr) *total_ptr = combined_map->assigned_score;
    
    return  0;
}
/***************************************/
/***************************************/

# if 0
    for ( uniq_map_ctr= 0;  uniq_map_ctr<list_uniq->map_max;  uniq_map_ctr++) {
	int i, j;
	current_map = list_uniq->map + uniq_map_ctr;
	printf ( "\n unique_map: %d  score: %8.3lf\n", uniq_map_ctr, current_map->assigned_score);
	for (i=0; i<current_map->x2y_size; i++) {
	     j = current_map->x2y[i];
	     printf ("  %3d  %3d    %8.3lf \n", i, j, current_map->sse_pair_score[i][j]);
	}
    }
    exit(1);
# endif
