/*
This source code is part of deconSTRUCT,
protein structure database search and backbone alignment application.
Written by Ivana Mihalek, with contributions from Mile Sikic.
Copyright (C) 2012 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or,
at your option, any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# include "struct.h"
# include "sys/time.h"
# ifdef OMP
#   include "omp.h"
# endif

	 //# include "gperftools/profiler.h"

# define  TOP_RMSD 200
# define  BAD_RMSD 10.0
# define JACKFRUIT 8
# define NUM_THREADS 8
# define CUTOFF_DNA 3.0

# define MAX_TRIPS  5000

typedef struct {
    int  member[3];
    int fingerprint; /* types of the three vectors and their hand */
} Triple;

# define TYPE1 1<<4
# define TYPE2 1<<3
# define TYPE3 1<<2
# define HAND  1<<1
# define UNDEF_HAND  1

int find_best_triples_exhaustive (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
				  double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
				  double **best_quat);
   
int find_best_triples_greedy(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
			     double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
			     double **best_quat);

int find_best_triples_exhaustive_parallel(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
        double **best_quat);

int find_best_triples_exhaustive_redux (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
				  double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
					double **best_quat);

int  sortTriplets(int ** best_triple_x_array, int ** best_triple_y_array,
		  double * best_rmsd_array, double ** best_quat_array, int top_rmsd);




/****************************************/
int complement_match (Representation* X_rep, Representation* Y_rep, List_of_maps *list){
	
    Penalty_parametrization penalty_params; /* for SW */
    double **x    = X_rep->full;
    int * x_type  = X_rep->full_type;
    int NX        = X_rep->N_full;
    double **y    = Y_rep->full;
    int * y_type  = Y_rep->full_type;
    int NY        = Y_rep->N_full;
    
    double F_effective = 0.0;
    double F_current;
    double q[4] = {0.0};
    double **x_rotated = NULL;
    double **tr_x_rotated = NULL;
    double **R;
    double z_scr = 0.0, *list_of_best_scores;
    double avg, avg_sq, stdev;
    double alpha = options.alpha;
    double rmsd, *best_rmsd;
    double **best_quat;
    double cutoff_rmsd = CUTOFF_DNA; /* <<<<<<<<<<<<<<<<< hardcoded */
    int *x_type_fudg, *y_type_fudg;
    int *anchor_x, *anchor_y, no_anchors;
    int no_top_rmsd = TOP_RMSD;
    int top_ctr;
    int **best_triple_x;
    int **best_triple_y;
    int retval; 
    int done = 0;
    int best_ctr;
    int i, j;
    int t;
    int smaller;
    int map_ctr  = 0;
    int no_maps_sorted;
    int * x2y, map_unstable;
    int current_map_id;
    //time_t  time_now, time_start;
    Map *current_map;
   
    int cull_by_dna (Representation * X_rep, int *set_of_directions_x,
		 Representation * Y_rep, int *set_of_directions_y,
		     int set_size, Map *map, double cutoff_rmsd);
    int distance_of_nearest_approach (Representation * X_rep, int *set_of_directions_x,
				      Representation * Y_rep, int *set_of_directions_y,
				      int set_size,  double * rmsd_ptr);
    int same_hand_triple (Representation * X_rep, int *set_of_directions_x,
			  Representation * Y_rep, int *set_of_directions_y, int set_size);
    
    int find_map (Penalty_parametrization * params, Representation *X_rep,  Representation *Y_rep,
		  double ** R, double alpha, double * F_effective,  Map *map, 
		  int *anchor_x, int * anchor_y, int anchor_size );
    int find_next_triple (double **X, double **Y, 
			  int *x_type, int *y_type, int NX, int NY,
			  int *x_triple, int *y_triple);
    int gradient_descent (int first_call, double alpha,
			  double **x, int *x_type, int NX,
			  double **y, int *y_type, int NY,
			  double *q_best, double *F_best_ptr) ;
    int monte_carlo (double alpha,
		 double **x, int * x_type, int NX,
		 double **y, int * y_type, int NY,
		 double  *q_best, double *F_best_ptr);
    int opt_quat (double ** x, int NX, int *set_of_directions_x,
		  double ** y, int NY, int *set_of_directions_y,
		  int set_size, double * q, double * rmsd);
    int qmap (double *x0, double *x1, double *y0, double *y1, double * quat);
    int store_sorted (List_of_maps * list,  double * best_score, int * new_map_id);
    
    smaller = (NX <= NY) ? NX : NY;
    no_top_rmsd = NX*NY/10; /* I'm not sure that this is the scale, but it works for now */
    if (no_top_rmsd < 100) no_top_rmsd = 100;

    
    /***********************/
    /* memory allocation   */
    /***********************/
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    if ( ! (x_rotated    = dmatrix (NX,3)) ) return 1;
    if ( ! (tr_x_rotated = dmatrix (NX,3)) ) return 1;
    if ( ! (best_quat    = dmatrix (no_top_rmsd,4)) ) return 1;
    if ( ! (best_rmsd    = emalloc (no_top_rmsd*sizeof(double))) ) return 1;
    if ( ! (best_triple_x    = intmatrix (no_top_rmsd,3)) ) return 1;
    if ( ! (best_triple_y    = intmatrix (no_top_rmsd,3)) ) return 1;
    if ( ! (list_of_best_scores = emalloc(list->no_maps_allocated*sizeof(double) )) ) return 1;
    if ( ! (x_type_fudg = emalloc(NX*sizeof(int) )) ) return 1;
    if ( ! (y_type_fudg = emalloc(NY*sizeof(int) )) ) return 1;
    if ( ! (anchor_x = emalloc(NX*sizeof(int) )) ) return 1;
    if ( ! (anchor_y = emalloc(NY*sizeof(int) )) ) return 1;

    /***********************/
    
    /***********************/
    /* expected quantities */
    /***********************/
    avg = avg_sq = stdev = 0.0;
    if (0) {
	if ( ! options.path[0] ) {
	    fprintf (stderr, "path to integral table must be given "
		     "in the cmd file (kwd \"path\") to calculate z-score for F.\n");
	    exit (1);
	}
	if (F_moments (x, x_type, NX, y, y_type, NY, alpha, &avg, &avg_sq, &stdev)) return 1;
    }
    /***********************/
    

    /***********************/
    /* initialization      */
    /***********************/
    best_ctr   = 0;
    penalty_params.gap_opening   = options.gap_open;
    penalty_params.gap_extension = options.gap_extend;
    penalty_params.endgap        = options.endgap;
    penalty_params.endgap_special_treatment = options.use_endgap;
    penalty_params.custom_gap_penalty_x = NULL;
    penalty_params.custom_gap_penalty_y = NULL;
    if ( ! (penalty_params.custom_gap_penalty_x = emalloc(NX*sizeof(double) )) ) return 1; 
    if ( ! (penalty_params.custom_gap_penalty_y = emalloc(NY*sizeof(double) )) ) return 1; 
    /***********************/


    /***************************************/
    /* find reasonable triples of SSEs     */
    /* that correspond in type             */
    /*  and can be mapped onto each other  */
    /***************************************/



    //ProfilerStart("profile.out") ; /* the output will got to the file called profile.out */
    if (options.exhaustive) {
	/*
	 * Exhaustive search for all triplets
	 */
	if (options.omp) {
# ifdef OMP
	    find_best_triples_exhaustive_parallel (X_rep, Y_rep, no_top_rmsd, best_rmsd, 
						   best_triple_x, best_triple_y, best_quat);
# else
	    infox ("To use omp, please recompile with -DOMP flag.\n",1);
# endif
	} else {
	    find_best_triples_exhaustive_redux (X_rep, Y_rep, no_top_rmsd, best_rmsd, 
	    				  best_triple_x, best_triple_y, best_quat);
	    //find_best_triples_exhaustive (X_rep, Y_rep, no_top_rmsd, best_rmsd, 
	    //best_triple_x, best_triple_y, best_quat);
	}
    } else {
    
	/*
	 * Greedy search - old algorithm
	 */
  	find_best_triples_greedy (X_rep, Y_rep, no_top_rmsd, best_rmsd,  
			      best_triple_x, best_triple_y, best_quat);
    }
    
    //ProfilerStop();
    current_map_id = 0;
    
    /*********************************************/
    /*   main loop                               */
    /*********************************************/
    no_maps_sorted = 0.0;
    map_ctr = 0;
    for (top_ctr=0; top_ctr<no_top_rmsd && done==0; top_ctr++) {

	if ( best_rmsd[top_ctr] > BAD_RMSD ) break;

	//printf (">> rmsd ok\n");
	
	quat_to_R (best_quat[top_ctr], R);
	rotate (x_rotated, NX, R, x);


	F_current = F( y, y_type, NY, x_rotated, x_type, NX, alpha);

	/* find map which uses the 2 triples as anchors */
	no_anchors = 3;
	
	current_map    = list->map+current_map_id;
	clear_map (current_map);
	find_map (&penalty_params, X_rep, Y_rep, R, alpha, &F_effective, current_map,
		   best_triple_x[top_ctr], best_triple_y[top_ctr], no_anchors);


	/* does this map still map the two triples we started with? */
	x2y = current_map->x2y;
	map_unstable  = 0;
	for (t=0; t<3; t++ ) {
	    if ( x2y[best_triple_x[top_ctr][t]] != best_triple_y[top_ctr][t] ) {
		map_unstable = 1;
	    }
	}
	if ( map_unstable) continue;

	//printf (">> map stable\n");

	/* do the mapped SSEs match in length? */
	if (options.use_length &&
	   current_map->avg_length_mismatch  > options.avg_length_mismatch_tol)  continue;
	
	//printf (">> passed length test\n");
	
	/* dna here is not DNA but "distance of nearest approach" */
	cull_by_dna ( X_rep, best_triple_x[top_ctr], 
		      Y_rep, best_triple_y[top_ctr],  3,  current_map, cutoff_rmsd );
	
	// why am I not dropping the thing if the rmsd id above cutoff?
	
	/* monte that optimizes the aligned vectors only */
	for (i=0; i<NX; i++) {
	     x_type_fudg[i] = JACKFRUIT;
	}
	for (j=0; j<NY; j++) {
	     y_type_fudg[j] = JACKFRUIT*2;
	}

	no_anchors = 0;

	for (i=0; i<NX; i++) {
	     j = current_map->x2y[i];
	     if (j < 0 ) continue;
	     x_type_fudg[i] = x_type[i];
	     y_type_fudg[j] = y_type[j];
	     anchor_x[no_anchors] = i;
	     anchor_y[no_anchors] = j;
	     no_anchors ++;
	}


	if ( opt_quat (x,  NX, anchor_x, y, NY, anchor_y, no_anchors, q, &rmsd)) continue;

	retval = monte_carlo (alpha, x, x_type_fudg, NX,  y,  y_type_fudg, NY, q, &F_current);



	if (retval) return retval;

	quat_to_R (q, R);
	/* store_sse_pair_score() is waste of time, but perhaps not critical */
	store_sse_pair_score (X_rep, Y_rep, R,  alpha, current_map);
	/* recalculate the assigned score*/
	map_assigned_score   (X_rep, current_map);


        /*   store the map that passed all the filters down to here*/
	current_map->F       = F_current;
	current_map->avg     = avg;
	current_map->avg_sq  = avg_sq;
	current_map->z_score = z_scr;
	memcpy ( current_map->q, q, 4*sizeof(double) );
	
	/************************/
	/* store sorted         */
	/************************/
	/* find the place for the new score and the new map */
    	store_sorted (list, list_of_best_scores, &current_map_id);
	
    }

    /**********************/
    /* garbage collection */
    gradient_descent (1, 0.0,  NULL, NULL, 0,
			       NULL, NULL, 0, NULL, NULL);
    free_dmatrix (R);
    free_dmatrix (x_rotated);
    free_dmatrix (tr_x_rotated);
    free_dmatrix (best_quat);
    free_imatrix (best_triple_x);
    free_imatrix (best_triple_y);
    free (list_of_best_scores);
    free (best_rmsd);
    free (x_type_fudg);
    free (y_type_fudg);
    free (anchor_x);
    free (anchor_y);
     
    if (penalty_params.custom_gap_penalty_x) free (penalty_params.custom_gap_penalty_x);
    if (penalty_params.custom_gap_penalty_y) free (penalty_params.custom_gap_penalty_y);
    /*********************/
   
    return 0;
}

/***************************************/
/***************************************/

int store_sorted (List_of_maps * list,  double * best_score, int * new_map_id) {

    Map *new_map =  list->map+ (*new_map_id); 
    int sorted_position, chunk, best_ctr = list->best_array_used;
    int map_max = list ->best_array_allocated;
    int *map_best = list->map_best;
    double new_score = -new_map->assigned_score;
    
    /* find the place in the list*/
    sorted_position = 0;
    while (sorted_position<best_ctr
	   && (new_score >= best_score[sorted_position])
	   && (sorted_position<map_max) )  sorted_position++;

    if ( sorted_position==map_max ) {/* we already have enough maps which are better */
	    
    } else { /* store the map  */
	
	
	/* store  the score and the map position*/
	/* move */
	chunk = best_ctr - sorted_position;
	if ( chunk ) {
	    memmove (best_score+sorted_position+1, best_score+sorted_position, chunk*sizeof(double)); 
	    memmove (  map_best+sorted_position+1,   map_best+sorted_position, chunk*sizeof(int)); 
	}
	best_score[sorted_position]   = new_score;
	map_best  [sorted_position]   = *new_map_id;

	/* what is the next available place in the map list? */
	if ( list->no_maps_used<map_max-1) {/* add to the first available place */
	    
	    list->no_maps_used++;
            list->best_array_used++;
	    *new_map_id = list->no_maps_used;
	    
	} else { /* add to the first available place */
	    int worst_map = map_best[map_max-1];
	    /* replace the worst map wiht the new  one */
	    *new_map_id = worst_map;
	}

    }

# if 0
    for ( sorted_position=0; sorted_position< list->no_maps_used; sorted_position++) {
	printf ( "  %3d   %8.3lf    0x%x    %8.3lf  %3d     0x%x     0x%x   \n",
		 sorted_position,  best_score[sorted_position],
		 list->map+ map_best[sorted_position], (list->map+ map_best[sorted_position])->assigned_score,
		 map_best[sorted_position], (list->map+ map_best[sorted_position])->cosine ,
		 (list->map+ map_best[sorted_position])->cosine[0] );
    }
    printf ("*****\n");
    infox ( "", 1);
# endif

    
    return 0;
}
/****************************************************************/
/****************************************************************/

int find_next_triple (double **X, double **Y, 
		    int *x_type, int *y_type,
		    int NX, int NY, int *x_triple, int *y_triple){

    int i, j, t, done = 0;
    int x_ctr;
    int y_ctr;
    
    /* next Y[y_ctr] cannot be numerically negative of X[x_ctr] */
    x_ctr = x_triple[0];
    y_ctr = y_triple[0];
    if ( x_type[x_ctr] != y_type[y_ctr] ) return 1;


    for (t=1; t<=2; t++ ) {
	x_ctr = x_triple[t-1];
	y_ctr = y_triple[t-1];

	if ( x_ctr == -1 || y_ctr==-1) break;
	
	x_triple[t] = -1;
	y_triple[t] = -1;
	done = 0;
	for (i=x_ctr+1; i<NX && !done; i++ ) {
	    for (j=y_ctr+1; j<NY && !done; j++ ) {
	
		if ( x_type[i] == y_type[j] ) {
		    x_triple[t] = i;
		    y_triple[t] = j;
		    done = 1;
		    
		}
	    }
	}
    }


    

    return 0;
}





/**********************************************************/
/**********************************************************/
/* check that the two are not mirror images - count on it being  a triple*/
int hand (Representation * X_rep,  int *set_of_directions_x) {
    
    double **x    = X_rep->full;
    double **x_cm = X_rep->cm;
    double cm_vector[3], avg_cm[3], cross[3], d, aux;
    int i, ray_a, ray_b,  ray_c, a, b, c;
    
 
 
    a = 0;
    b = 1;
    c = 2;
    
    /*****************************/
    /*****************************/
    ray_a = set_of_directions_x[a];
    ray_b = set_of_directions_x[b];
    ray_c = set_of_directions_x[c];
    /* I better not use the cross prod here: a small diff
       int he angle making them pointing toward each other or
       away from each other changes the direction of cross prod;
       rather, use one vector, and distance between the cm's as the other */
    for (i=0; i<3; i++ ) {
	cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
    }
    if ( ! normalized_cross (x[ray_a], cm_vector, cross, &aux) ) {
    
	/* note I am making another cm vector here */
	for (i=0; i<3; i++ ) {
	    avg_cm[i] = (x_cm[ray_b][i] + x_cm[ray_a][i])/2;
	    cm_vector[i] = x_cm[ray_c][i] - avg_cm[i];
	}
	unnorm_dot (cm_vector, cross, &d);

     
	if ( d > 0 ) {
	    return 0;
	} else {
	    return HAND;
	}
	
    } else {
	/* the cross product cannot be calculated,
	   presumably bcs the vectors are (anti) parallel */
	return UNDEF_HAND;
    }
    
}

/**********************************************************/
/**********************************************************/
/* check that the two are not mirror images - count on it being  a triple*/
int same_hand_triple (Representation * X_rep,  int *set_of_directions_x,
	       Representation * Y_rep, int *set_of_directions_y, int set_size) {
    
    double **x    = X_rep->full;
    double **x_cm = X_rep->cm;
    double **y    = Y_rep->full;
    double **y_cm = Y_rep->cm;
    double cm_vector[3], avg_cm[3], cross[3],  dx, dy, aux;
    int i, ray_a, ray_b,  ray_c, a, b, c;
    
    if (set_size !=3 ) return 0;

 
    a = 0;
    b = 1;
    c = 2;
    
    /*****************************/
    /*****************************/
    ray_a = set_of_directions_x[a];
    ray_b = set_of_directions_x[b];
    ray_c = set_of_directions_x[c];
    
    /* I better not use the cross prod here: a small diff
       int he angle makeing them pointing toward each other or
       away from each other changes the direction of cross prod;
       rather, use one vector, and distance between the cm's as the other */
    for (i=0; i<3; i++ ) {
	cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
    }


    /* note that the cross product here is between the direction
       vector of SSE labeled a, and the vector spanning the geeometric
       centers of a and b - this cross prod can be undefined only if
       the two are stacked */
    if ( ! normalized_cross (x[ray_a], cm_vector, cross, &aux)) {
	/* the function returns 0 on success, i.e. retval==0 means
	   the cross product can be calculated ok. */ 
	
	/* note I am making another cm vector here */
	for (i=0; i<3; i++ ) {
	    avg_cm[i] = (x_cm[ray_b][i] + x_cm[ray_a][i])/2;
	    cm_vector[i] = x_cm[ray_c][i] - avg_cm[i];
	}
	unnorm_dot (cm_vector, cross, &dx);

	/*****************************/
	/*****************************/
	ray_a = set_of_directions_y[a];
	ray_b = set_of_directions_y[b];
	ray_c = set_of_directions_y[c];
	for (i=0; i<3; i++ ) {
	    cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
	}
	if ( !normalized_cross (y[ray_a], cm_vector, cross, &aux)) {

	    /* if  I can calculate the cross product for one triplet
	       but not for the other, I declare it a different hand (this will
	       be dropped from further consideration) */
	    return 0;

	}
	/* note I am makning another cm vector here */
	for (i=0; i<3; i++ ) {
	    avg_cm[i] = (y_cm[ray_b][i] + y_cm[ray_a][i])/2;
	    cm_vector[i] = y_cm[ray_c][i] - avg_cm[i];
	}
	/*note: unnorm_dot thinks it is getting unit vectors,
	  and evrything that is >1 will be "rounded" to 1
	  (similarly for -1) - it doesn't do the normalization itself*/
	unnorm_dot (cm_vector, cross, &dy);

     
	if ( dx*dy < 0 ) {
	    return 0;
	} else {
 	    return 1;   /* this isn't err value - the handedness is the same */
	}

	
    } else {

	/* if the vectors linking cms are parallel, the hand is the same */
	double cm_vector_2[3], dot;
	
	ray_a = set_of_directions_y[a];
	ray_b = set_of_directions_y[b];
	ray_c = set_of_directions_y[c];
	for (i=0; i<3; i++ ) {
	    cm_vector_2[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
	}

	unnorm_dot(cm_vector, cm_vector_2, &dot);

	if ( dot > 0) {
	    return 1;
	} else {
	    return 0;
	}

	/* (hey, at least I am checking the retvals of my functions) */
	
    }



    return 0; /* better safe than sorry */
}

/**********************************************************/

int distance_of_nearest_approach ( Representation * X_rep,  int *set_of_directions_x,
				   Representation * Y_rep, int *set_of_directions_y,
				   int set_size,  double * rmsd_ptr) {
    
    double **x    = X_rep->full;
    double **x_cm = X_rep->cm;
    double **y    = Y_rep->full;
    double **y_cm = Y_rep->cm;
    double cm_vector[3], cross[3],  distance_x, distance_y;
    double aux, rmsd;
    int i, ray_a, ray_b, a, b, prev_next, norm;
    
    if (set_size <=1 ) return 1;

    rmsd = 0.0;
    norm = 0;
    
    /* the rmsd for the remaining vectors is ... */
    for (a=0; a<set_size; a++) {

	for ( prev_next = -1;  prev_next <= +1; prev_next +=2 ) {
	    b = (set_size+a+prev_next)%set_size;
	
	    ray_a = set_of_directions_x[a];
	    ray_b = set_of_directions_x[b];
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions x */
	    for (i=0; i<3; i++ ) {
		cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
	    }
	    
	    if ( normalized_cross (x[ray_a], cm_vector, cross, &distance_x) ) {
		
		/* x[ray_a], cm_vector  are (anti) parallel, corss_product undefined */
		/* we then expect the same kind of behavior fromt he triple in the
		   other strcuture*/
		distance_x = 0;
		for (i=0; i<3; i++ ) {
		    distance_x += cm_vector[i]*cm_vector[i];
		}
		distance_x = sqrt(distance_x);
		
		ray_a = set_of_directions_y[a];
		ray_b = set_of_directions_y[b];
	    
		/* distance of nearest approach of ray b
		   to the cm of a, in the set of directions y */
		for (i=0; i<3; i++ ) {
		    cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
		}
	    
		distance_y = 0;
		for (i=0; i<3; i++ ) {
		    distance_y += cm_vector[i]*cm_vector[i];
		}
		distance_y = sqrt(distance_y);

	    } else {

		ray_a = set_of_directions_y[a];
		ray_b = set_of_directions_y[b];
	    
		/* distance of nearest approach of ray b
		   to the cm of a, in the set of directions y */
		for (i=0; i<3; i++ ) {
		    cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
		}
	    
		if ( normalized_cross (cm_vector, y[ray_a], cross, &distance_y) ) {
		    /* one triplet the vectors are parallel, in the other one they are not ...*/
		    /* but that's numerical ...*/
		    /* I'll just add some penalty here: */
		    distance_y = distance_x + CUTOFF_DNA*1.1;
		}
	    }
		
	    /* distance_x and distance_y shold be about the same, if these are two equivalent triples */
	    aux   = distance_x-distance_y;
	    rmsd += aux*aux;
	    norm ++;
# if 0
	    printf ("%d, %d   x:  %2d  %2d  y:  %2d  %2d  \n", a, b,
		    set_of_directions_x[a], set_of_directions_x[b], 
		    set_of_directions_y[a], set_of_directions_y[b]); 
	    printf (" distance x:  %8.3lf  distance y:  %8.3lf   difference:   %8.3lf \n",
		    distance_x,  distance_y, fabs (distance_x-distance_y));
# endif
	}
	
    }

    rmsd /= norm;
    rmsd  = sqrt(rmsd);
    *rmsd_ptr = rmsd;
    
    return 0;
}

/**********************************************************/
/**********************************************************/
int cull_by_dna (Representation * X_rep, int *set_of_directions_x,
		 Representation * Y_rep, int *set_of_directions_y,
		 int set_size, Map *map, double cutoff_rmsd) {
    
    double **x    = X_rep->full;
    double **x_cm = X_rep->cm;
    double **y    = Y_rep->full;
    double **y_cm = Y_rep->cm;
    double cm_vector[3], cross[3],  distance_x, distance_y;
    double aux, rmsd;
    int NX=X_rep->N_full;
    int i, j, ray_a, ray_b, a, b, prev_next;
    int in_triple, norm;
    
    if (set_size < 1 ) return 1;


    for (i=0; i<NX; i++) {
	j = map->x2y[i];
	if ( j<0) continue;

	/* check that i is not in the triple here */
	in_triple = 0;
	for (a=0; a<set_size; a++) {
	    if ( i==set_of_directions_x[a]) {
		in_triple =1;
		break;
	    }
	}
	if (in_triple) continue;



	rmsd = 0.0;
	norm = 0;
    
	/* the rmsd for the remaining vectors is ... */
	for (a=0; a<set_size; a++) {

	    for ( prev_next = -1;  prev_next <= +1; prev_next +=2 ) {
		b = (set_size+a+prev_next)%set_size;
	
		ray_a = set_of_directions_x[a];
		ray_b = set_of_directions_x[b];
		/* distance of nearest approach of ray b
		   to the cm of a, in the set of directions x */
		for (i=0; i<3; i++ ) {
		    cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
		}
	    
		if ( normalized_cross (x[ray_a], cm_vector, cross, &distance_x) ) {
		
		    /* x[ray_a], cm_vector  are (anti) parallel, corss_product undefined */
		    /* we then expect the same kind of behavior fromt he triple in the
		       other strcuture*/
		    distance_x = 0;
		    for (i=0; i<3; i++ ) {
			distance_x += cm_vector[i]*cm_vector[i];
		    }
		    distance_x = sqrt(distance_x);
		
		    ray_a = set_of_directions_y[a];
		    ray_b = set_of_directions_y[b];
	    
		    /* distance of nearest approach of ray b
		       to the cm of a, in the set of directions y */
		    for (i=0; i<3; i++ ) {
			cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
		    }
	    
		    distance_y = 0;
		    for (i=0; i<3; i++ ) {
			distance_y += cm_vector[i]*cm_vector[i];
		    }
		    distance_y = sqrt(distance_y);

		} else {

		    ray_a = set_of_directions_y[a];
		    ray_b = set_of_directions_y[b];
	    
		    /* distance of nearest approach of ray b
		       to the cm of a, in the set of directions y */
		    for (i=0; i<3; i++ ) {
			cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
		    }
	    
		    if ( normalized_cross (cm_vector, y[ray_a], cross, &distance_y) ) {
			/* one triplet the vectors are parallel, in the other one they are not ...*/
			/* but that's numerical ...*/
			/* I'll just add some penalty here: */
			distance_y = distance_x + CUTOFF_DNA*1.1;
		    }
		}
		
		/* distance_x and distance_y shold be about the same, if these are two equivalent triples */
		aux   = distance_x-distance_y;
		rmsd += aux*aux;
		norm ++;
# if 0
		printf ("%d, %d   x:  %2d  %2d  y:  %2d  %2d  \n", a, b,
			set_of_directions_x[a], set_of_directions_x[b], 
			set_of_directions_y[a], set_of_directions_y[b]); 
		printf (" distance x:  %8.3lf  distance y:  %8.3lf   difference:   %8.3lf \n",
			distance_x,  distance_y, fabs (distance_x-distance_y));
# endif
	    }
	
	}

	rmsd /= norm;
	rmsd  = sqrt(rmsd);




        /* if rmsd is bigger than the cutoff, lose this from the mapping */
	if ( rmsd > cutoff_rmsd) {
	    map->x2y[i] = -1;
	    map->y2x[j] = -1;
	}
    }
    
     
    return 0;
}

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
	exit (1);
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
    
    if ( !( A=dmatrix(4,4) ) ) exit (1);
    memcpy (A[0], ATA_sum[0], 4*4*sizeof(double));


   /* note how we pass the matrix: */
    dsyev_ ( &jobz, &uplo,  &n, A[0], &lda, w, work, &lwork, &info);
    if (  ! info) {
	*rmsd = sqrt (w[0]);
	for (i=0; i<4; i++ ) q[i] = A[0][i];
	if (0) {
	    /* w contains the eigenvalues */
	    printf ("\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", w[i]);
	    printf ("\nrmsd: %8.3lf \n", *rmsd);
	    printf ("quat:\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", q[i]);
	    printf ("\n");
	    /* printf (" opt lwork: %d\n", (int) work[0]); */
	}
    } else {
	fprintf (stderr, "Error in dsyev().\n");
	exit (1);
    }
    
   
    
    free_dmatrix(A);

    return 0;
    
}

/************************************/

/************************************/


/***************************************/
/* find reasonable triples of SSEs      */
/* that correspond in type             */
/*  and can be mapped onto each other  */
/***************************************/

/**
 * 
 * @param X_rep
 * @param Y_rep
 * @param no_top_rmsd
 * @param best_rmsd
 * @param best_triple_x
 * @param best_triple_y
 * @param best_quat
 * @return
 *
 */


int find_best_triples_exhaustive_redux (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
				  double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
				  double **best_quat) {

    int top_ctr, i_x, j_x, k_x, i_y, j_y, k_y;
    int xtrip_ct, ytrip_ct;
    int no_xtrips, no_ytrips, no_trip_pairs_to_compare;
    int y_list_full, x_list_full;
    int x_enumeration_done, y_enumeration_done;
    int panic_ctr;
    int chunk, same, i;
    int NX = X_rep->N_full;
    int NY = Y_rep->N_full;
    int * x_type = X_rep->full_type;
    int * y_type = Y_rep->full_type;
    Triple *x_triple, *y_triple;
 
    double **x = X_rep->full;    
    double **y = Y_rep->full;

    double epsilon = 0.05;
    double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
    double rmsd;
    double q_init[4] = {0.0};
    double threshold_dist = options.threshold_distance;
     
    double ** cmx = X_rep->cm;
    double ** cmy = Y_rep->cm;

    if (options.verbose) printf ("exhaustive search \n");

    /***************************************/
    /* find reasonable triples of SSEs      */
    /* that correspond in type             */
    /*  and can be mapped onto each other  */
    /***************************************/
    for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
        best_rmsd[top_ctr] = BAD_RMSD + 1;
        best_triple_x[top_ctr][0] = -1;
    }

    if ( ! (x_triple = emalloc (MAX_TRIPS*sizeof(Triple) ) ) )  return 1;
    if ( ! (y_triple = emalloc (MAX_TRIPS*sizeof(Triple) ) ) )  return 1;
    

    no_trip_pairs_to_compare = 0;
    x_enumeration_done = 0;
    i_x = 0;
    j_x = 1;
    k_x = 1;

    panic_ctr = 0;	
    while ( ! x_enumeration_done ) {
	
	x_list_full = 0;
	xtrip_ct = 0;

	while ( ! x_list_full  && !x_enumeration_done) {

	    k_x ++;
	    if ( k_x == NX ) {
		    
		if ( j_x <  NX-2 ) {
		    j_x++;
		    k_x = j_x+1;
			
		} else {
			
		    if ( i_x <  NX-3) {
			i_x++;
			j_x = i_x+1;
			k_x = j_x+1;
			    
		    } else {
			x_enumeration_done = 1;
			break; /* from filling the x list */
		    }			
		}
	    }
	    
	    if (two_point_distance(cmx[i_x],cmx[j_x]) > threshold_dist) continue;
	    if (two_point_distance(cmx[i_x],cmx[k_x]) > threshold_dist) continue;  
	    if (two_point_distance(cmx[j_x],cmx[k_x]) > threshold_dist) continue;  
	    		
	    x_triple[xtrip_ct].fingerprint = 0;
	    if (x_type[i_x] == HELIX) x_triple[xtrip_ct].fingerprint |= TYPE1;
	    if (x_type[j_x] == HELIX) x_triple[xtrip_ct].fingerprint |= TYPE2;
	    if (x_type[k_x] == HELIX) x_triple[xtrip_ct].fingerprint |= TYPE3;
	    

	    x_triple[xtrip_ct].member[0] = i_x;
	    x_triple[xtrip_ct].member[1] = j_x;
	    x_triple[xtrip_ct].member[2] = k_x;

	    x_triple[xtrip_ct].fingerprint |=  hand(X_rep, x_triple[xtrip_ct].member);
		
	    xtrip_ct++;
	    x_list_full = ( xtrip_ct == MAX_TRIPS);

	    
	} /* end filling the x trips list */


	
	no_xtrips = xtrip_ct;

	 
	y_enumeration_done = 0;
	i_y = 0;
	j_y = 1;
	k_y = 1;
	while ( !y_enumeration_done) {


	    y_list_full = 0;
	    ytrip_ct = 0;
	    
	    while ( ! y_list_full  && !y_enumeration_done) {

		k_y ++;
		if ( k_y == NY ) {
		    
		    if ( j_y <  NY-2 ) {
			j_y++;
			k_y = j_y+1;
			
		    } else {
			
			if ( i_y <  NY-3) {
			    i_y++;
			    j_y = i_y+1;
			    k_y = j_y+1;
			    
			} else {
			    y_enumeration_done = 1;
			    break; /* from filling the y list */
			}			
		    }
		}
		
		if (two_point_distance(cmy[i_y],cmy[j_y]) > threshold_dist) continue;
		if (two_point_distance(cmy[i_y],cmy[k_y]) > threshold_dist) continue;  
		if (two_point_distance(cmy[j_y],cmy[k_y]) > threshold_dist) continue;
		
		y_triple[ytrip_ct].fingerprint = 0;
		if (y_type[i_y] == HELIX) y_triple[ytrip_ct].fingerprint |= TYPE1;
		if (y_type[j_y] == HELIX) y_triple[ytrip_ct].fingerprint |= TYPE2;
		if (y_type[k_y] == HELIX) y_triple[ytrip_ct].fingerprint |= TYPE3;
		
	    
		
		y_triple[ytrip_ct].member[0] = i_y;
		y_triple[ytrip_ct].member[1] = j_y;
		y_triple[ytrip_ct].member[2] = k_y;

		if ( hand(Y_rep, y_triple[ytrip_ct].member ) ) y_triple[ytrip_ct].fingerprint |= HAND;
		
		ytrip_ct++;
		y_list_full = ( ytrip_ct == MAX_TRIPS);
	    
	    } /* end filling the y trips list */

	    
	    no_ytrips = ytrip_ct;
	    panic_ctr ++;
	    if ( panic_ctr == 1000 ) {
		fprintf (stderr, "Exit by panic ctr in %s:%d.\n", __FILE__, __LINE__);
		exit (1);
	    }
	
	    //printf ("\t no trips in x: %d\n",  no_xtrips);
	    //printf ("\t no trips in y: %d\n\n",  no_ytrips);
	    for (xtrip_ct=0; xtrip_ct<no_xtrips; xtrip_ct++) {
		for (ytrip_ct=0; ytrip_ct<no_ytrips; ytrip_ct++) {

		    no_trip_pairs_to_compare ++;
		    
		    /* filters :*/
		    if ( x_triple[xtrip_ct].fingerprint ^ y_triple[ytrip_ct].fingerprint) continue;

		    if (distance_of_nearest_approach(X_rep, x_triple[xtrip_ct].member,
						     Y_rep, y_triple[ytrip_ct].member, 3, &rmsd)) continue;
		    if ( rmsd > cutoff_rmsd) continue;
	    
		    if (opt_quat(x, NX, x_triple[xtrip_ct].member, y, NY, y_triple[ytrip_ct].member, 3, q_init, &rmsd)) continue;

		    same = 0;
		    for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
			    
			if (best_rmsd[top_ctr] > BAD_RMSD) break;
			same = 1;
			for (i=0; i<4; i++) {
			    if ( fabs(q_init[i] - best_quat[top_ctr][i]) > epsilon) {
				same = 0;
				break;
			    }
			}
			if ( same )  break;
			
				
		    }
		    if (same) continue;
			    

		    /* store the survivors */
		    for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {

			if (rmsd <= best_rmsd[top_ctr]) {

			    chunk = no_top_rmsd - top_ctr - 1;

			    if (chunk) {
				memmove(best_rmsd + top_ctr + 1, best_rmsd + top_ctr, chunk * sizeof (double));
				memmove(best_quat[top_ctr + 1],
					best_quat[top_ctr], chunk * 4 * sizeof (double));
				memmove(best_triple_x[top_ctr + 1],
					best_triple_x[top_ctr], chunk * 3 * sizeof (int));
				memmove(best_triple_y[top_ctr + 1],
					best_triple_y[top_ctr], chunk * 3 * sizeof (int));
			    }
			    best_rmsd[top_ctr] = rmsd;
			    memcpy(best_quat[top_ctr], q_init, 4 * sizeof (double));
			    memcpy(best_triple_x[top_ctr], x_triple[xtrip_ct].member, 3 * sizeof (int));
			    memcpy(best_triple_y[top_ctr], y_triple[ytrip_ct].member, 3 * sizeof (int));
                                    
			    break;

			}
		    }

		}
	    }

	    
	} /* y enumeration loop */
    } /* x enumeration loop  */
    

# if 0
    printf ("no trips in x: %d\n",  no_xtrips);
    printf ("no trips in y: %d\n",  no_ytrips);
    printf ("no trips to compare : %d\n",  no_trip_pairs_to_compare);
    
   for (top_ctr=0; top_ctr<no_top_rmsd; top_ctr++) {

	if ( best_rmsd[top_ctr] > BAD_RMSD ) break;
	printf (" %8.3lf      %3d  %3d  %3d     %3d  %3d  %3d \n",   best_rmsd[top_ctr],
		best_triple_x[top_ctr][0], best_triple_x[top_ctr][1], best_triple_x[top_ctr][2], 
		best_triple_y[top_ctr][0], best_triple_y[top_ctr][1], best_triple_y[top_ctr][2]);
   }

   // exit (1);
# endif
    


    free (x_triple);
    free (y_triple);
    
    return 0;

}


/**
 * 
 * @param X_rep
 * @param Y_rep
 * @param no_top_rmsd
 * @param best_rmsd
 * @param best_triple_x
 * @param best_triple_y
 * @param best_quat
 * @return 
 */

int find_best_triples_exhaustive (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
				  double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
				  double **best_quat) {

    int top_ctr, i, j, k, l, m, n;
    double **x = X_rep->full;
    int * x_type = X_rep->full_type;
    int NX = X_rep->N_full;
    double **y = Y_rep->full;
    int * y_type = Y_rep->full_type;
    int NY = Y_rep->N_full;
    int x_triple[3], y_triple[3];
    int chunk, same;
    double epsilon = 0.05;
    double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
    double rmsd;
    double q_init[4] = {0.0};
    double ** cmx = X_rep->cm;
    double ** cmy = Y_rep->cm;
    double threshold_dist = options.threshold_distance;

    if (options.verbose) printf ("exhaustive search \n");

    /***************************************/
    /* find reasonable triples of SSEs     */
    /* that correspond in type             */
    /*  and can be mapped onto each other  */
    /***************************************/
    for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
        best_rmsd[top_ctr] = BAD_RMSD + 1;
        best_triple_x[top_ctr][0] = -1;
    }

    /*
     * Exhaustive search through a 6D space - ugly code
     */

    for (i = 0; i < NX - 2; ++i) {
        x_triple[0] = i;

        for (j = 0; j < NY - 2; ++j) {
            if (x_type[i] != y_type[j]) continue;
            y_triple[0] = j;

           for (k = i + 1; k < NX -1 ; ++k) {
                if (two_point_distance(cmx[i],cmx[k]) > threshold_dist) continue;
                x_triple[1] = k;

               for (l = j + 1; l < NY -1 ; ++l) {
                    if (x_type[k] != y_type[l]) continue;
                    if (two_point_distance(cmy[j],cmy[l]) > threshold_dist) continue;
                    y_triple[1] = l;
		    
                    for (m = k + 1; m < NX; ++m) {
                        if (two_point_distance(cmx[i],cmx[m]) > threshold_dist) continue;
                        if (two_point_distance(cmx[k],cmx[m]) > threshold_dist) continue;
                        x_triple[2] = m;

                       for (n = l + 1; n < NY; ++n) {
                            if (two_point_distance(cmy[j],cmy[n]) > threshold_dist) continue;
                            if (two_point_distance(cmy[l],cmy[n]) > threshold_dist) continue;
                            if (x_type[m] != y_type[n]) continue;
                            y_triple[2] = n;


			    if (!same_hand_triple(X_rep, x_triple, Y_rep, y_triple, 3)) continue;

			    if (distance_of_nearest_approach(X_rep, x_triple,
                                    Y_rep, y_triple, 3, &rmsd)) continue;
                            
                            if (rmsd > cutoff_rmsd) continue;

 	
			
			    if (opt_quat(x, NX, x_triple, y, NY, y_triple, 3, q_init, &rmsd)) continue;

			    /* do we already have a similar quat ?*/
			    
			    for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
			    
				if (best_rmsd[top_ctr] > BAD_RMSD) break;
				same = 1;
				for (i=0; i<4; i++) {
				    if ( fabs(q_init[i] - best_quat[top_ctr][i]) > epsilon) {
					same = 0;
					break;
				    }
				}
				if ( same ) {
				    printf ("same  %8.3f  %8.3f      %8.3f  %8.3f    %8.3f  %8.3f     %8.3f  %8.3f    \n",
					    q_init[0],  best_quat[top_ctr][0],
					    q_init[1],  best_quat[top_ctr][1],
					    q_init[2],  best_quat[top_ctr][2],
					    q_init[3],  best_quat[top_ctr][3]
					);
				    printf ("   %3d %3d %3d     %3d %3d %3d    %8.3lf  \n",
					    best_triple_x[top_ctr][0],  best_triple_x[top_ctr][1],  best_triple_x[top_ctr][2],
					    best_triple_y[top_ctr][0],  best_triple_y[top_ctr][1],  best_triple_y[top_ctr][2],
					    best_rmsd[top_ctr] );
				    
				    printf ("   %3d %3d %3d     %3d %3d %3d   \n\n",
					    x_triple[0],  x_triple[1],  x_triple[2],
					    y_triple[0],  y_triple[1],  y_triple[2]);
				    
				    exit (1);
				}
				
			    }
		    
			    if ( x_triple[0] == 10 &&
				  x_triple[1] == 11 &&
				  x_triple[2] == 12 &&
				  y_triple[0] ==  3 &&
				  y_triple[1] ==  4 &&
				 y_triple[2] == 5 ) {
				printf ( "blah  %8.3f\n", rmsd);
				
				
			    }
			    for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {

                                if (rmsd <= best_rmsd[top_ctr]) {

				    chunk = no_top_rmsd - top_ctr - 1;

                                    if (chunk) {
                                        memmove(best_rmsd + top_ctr + 1, best_rmsd + top_ctr, chunk * sizeof (double));
                                        memmove(best_quat[top_ctr + 1],
                                                best_quat[top_ctr], chunk * 4 * sizeof (double));
                                        memmove(best_triple_x[top_ctr + 1],
                                                best_triple_x[top_ctr], chunk * 3 * sizeof (int));
                                        memmove(best_triple_y[top_ctr + 1],
                                                best_triple_y[top_ctr], chunk * 3 * sizeof (int));
                                    }
                                    best_rmsd[top_ctr] = rmsd;
                                    memcpy(best_quat[top_ctr], q_init, 4 * sizeof (double));
                                    memcpy(best_triple_x[top_ctr], x_triple, 3 * sizeof (int));
                                    memcpy(best_triple_y[top_ctr], y_triple, 3 * sizeof (int));
                                    
                                    break;

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    return 0;

}


/**
 * A parallel algorithm for exhaustive search of all triplets combinations using OpenMP
 * @param X_rep
 * @param Y_rep
 * @param no_top_rmsd
 * @param best_rmsd
 * @param best_triple_x
 * @param best_triple_y
 * @param best_quat
 * @return 
 */


int find_best_triples_exhaustive_parallel(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
        double **best_quat) {
# ifdef OMP
    // initialization of global array of values
    double ** best_quat_array = dmatrix(no_top_rmsd * NUM_THREADS, 4);
    int ** best_triple_x_array = intmatrix(no_top_rmsd * NUM_THREADS, 3);
    int ** best_triple_y_array = intmatrix(no_top_rmsd * NUM_THREADS, 3);
    double * best_rmsd_array = (double *) malloc(no_top_rmsd * NUM_THREADS * sizeof (double));
    
    int cnt;
    
    for (cnt = 0; cnt < NUM_THREADS*no_top_rmsd; ++cnt) {
        best_rmsd_array[cnt] = BAD_RMSD + 1;
        best_triple_x_array[cnt][0] = -1;
    }
        


   omp_set_num_threads(NUM_THREADS);

#pragma omp parallel 
    {

        int top_ctr, i, j, k, l, n;
        int myid = omp_get_thread_num();

        double ** best_quat_local = dmatrix(no_top_rmsd, 4);
        int ** best_triple_x_local = intmatrix(no_top_rmsd, 3);
        int ** best_triple_y_local = intmatrix(no_top_rmsd, 3);
        double * best_rmsd_local = (double *) malloc(no_top_rmsd * sizeof (double));
        double **x = X_rep->full; // no change
        int * x_type = X_rep->full_type; // no change
        int NX = X_rep->N_full; // no change
        double **y = Y_rep->full;
        int * y_type = Y_rep->full_type;
        int NY = Y_rep->N_full;
        int x_triple[3], y_triple[3];
        int chunk;
        double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
        double rmsd; // 
        double q_init[4] = {0.0}; // no change
        double ** cmx = X_rep->cm; // no change
        double ** cmy = Y_rep->cm; // no change
        double threshold_dist = options.threshold_distance;

        /***************************************/
        /* find reasonable triples of SSEs      */
        /* that correspond in type             */
        /*  and can be mapped onto each other  */
        /***************************************/
        for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
            best_rmsd_local[top_ctr] = BAD_RMSD + 1;
            best_triple_x_local[top_ctr][0] = -1;
        }

        /*
         * Exhaustive search through a 6D space - ugly code
         * Parallelization 
         */

#pragma omp for       
        for (i = 0; i < NX - 2; ++i) {
            int m;         
            
            x_triple[0] = i;
            for (j = 0; j < NY - 2; ++j) {
                if (x_type[i] != y_type[j]) continue;
                y_triple[0] = j;
                for (k = i + 1; k < NX - 1; ++k) {
                   if (two_point_distance(cmx[i], cmx[k]) > threshold_dist) continue;
                    x_triple[1] = k;
                    for (l = j + 1; l < NY - 1; ++l) {
                        if (x_type[k] != y_type[l]) continue;
                        if (two_point_distance(cmy[j], cmy[l]) > threshold_dist) continue;
                        y_triple[1] = l;
                        for (m = k + 1; m < NX; ++m) {
                            if (two_point_distance(cmx[i], cmx[m]) > threshold_dist) continue;
                            if (two_point_distance(cmx[k], cmx[m]) > threshold_dist) continue;
                            x_triple[2] = m;
                            for (n = l + 1; n < NY; ++n) {
                                if (two_point_distance(cmy[j], cmy[n]) > threshold_dist) continue;
                                if (two_point_distance(cmy[l], cmy[n]) > threshold_dist) continue;
                                if (x_type[m] != y_type[n]) continue;
                                y_triple[2] = n;

                                if (!same_hand_triple(X_rep, x_triple, Y_rep, y_triple, 3)) continue;

                                if (distance_of_nearest_approach(X_rep, x_triple,
                                        Y_rep, y_triple, 3, &rmsd)) continue;

                                if (rmsd > cutoff_rmsd) continue;
                                

                                if (opt_quat(x, NX, x_triple, y, NY, y_triple, 3, q_init, &rmsd)) continue;
                                for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
                                    // insertion of a new values in arrays keeping arrays sorted
                                    
                                    if (rmsd <= best_rmsd_local[top_ctr]) {
                                        chunk = no_top_rmsd - top_ctr - 1;

                                        if (chunk) {
                                            memmove(best_rmsd_local + top_ctr + 1,
						    best_rmsd_local + top_ctr, chunk * sizeof (double));
                                            memmove(best_quat_local[top_ctr + 1],
						    best_quat_local[top_ctr], chunk * 4 * sizeof (double));
                                            memmove(best_triple_x_local[top_ctr + 1],
                                                    best_triple_x_local[top_ctr], chunk * 3 * sizeof (int));
                                            memmove(best_triple_y_local[top_ctr + 1],
                                                    best_triple_y_local[top_ctr], chunk * 3 * sizeof (int));
                                        }
                                        best_rmsd_local[top_ctr] = rmsd;
                                        memcpy(best_quat_local[top_ctr], q_init, 4 * sizeof (double));
                                        memcpy(best_triple_x_local[top_ctr], x_triple, 3 * sizeof (int));
                                        memcpy(best_triple_y_local[top_ctr], y_triple, 3 * sizeof (int));

                                        break;

                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            // each thread copies values to global arrays in accordance with its thread id
            memcpy(*(best_quat_array + myid*no_top_rmsd), *(best_quat_local), no_top_rmsd * 4 * sizeof(double));
            memcpy(*(best_triple_y_array + myid*no_top_rmsd), *(best_triple_y_local), no_top_rmsd * 3 * sizeof(int));
            memcpy(*(best_triple_x_array + myid*no_top_rmsd), *(best_triple_x_local), no_top_rmsd * 3 * sizeof(int));
            memcpy(best_rmsd_array + myid*no_top_rmsd, best_rmsd_local, no_top_rmsd * sizeof(double));
            
            
        }
        

        free_dmatrix(best_quat_local);
        free_imatrix(best_triple_x_local);
        free_imatrix(best_triple_y_local);
        free(best_rmsd_local);  
        
        // parallel sort of elements of arrays 
        sortTriplets(best_triple_x_array, best_triple_y_array, best_rmsd_array, best_quat_array, no_top_rmsd);
    }
     
   // 
    memcpy(*best_quat, *best_quat_array, no_top_rmsd * 4 * sizeof(double));
    memcpy(*best_triple_y, *best_triple_y_array, no_top_rmsd * 3 * sizeof(int));
    memcpy(*best_triple_x, *best_triple_x_array, no_top_rmsd * 3 * sizeof(int));
    memcpy(best_rmsd, best_rmsd_array, no_top_rmsd * sizeof(double));
    

    free_dmatrix(best_quat_array);
    free_imatrix(best_triple_x_array);
    free_imatrix(best_triple_y_array);
    free(best_rmsd_array);

# endif
    
    return 0;

}


/**
 * A parallel algorithm for sorting triplets using OpenMP
 * @param best_triple_x_array
 * @param best_triple_y_array
 * @param best_rmsd_array
 * @param best_quat_array
 * @return 
 */
 
int  sortTriplets(int ** best_triple_x_array, int ** best_triple_y_array,
		  double * best_rmsd_array, double ** best_quat_array, int top_rmsd){

# ifdef OMP
    
     int stride, j, k, chunk;
     int myid = omp_get_thread_num();
 
     for (stride = NUM_THREADS/2; stride > 0; stride /= 2) {
         #pragma omp barrier
         j = 0;
         if (myid < stride){
             for (k = 0; k < top_rmsd; ++k){                
                 while (j < top_rmsd) {
                     if (best_rmsd_array[myid*top_rmsd + stride*top_rmsd  + k]
			 < best_rmsd_array[myid*top_rmsd + j]) {
                         
                         chunk = top_rmsd - j - 1;

                         if (chunk) {
                             memmove(best_rmsd_array + myid*top_rmsd + j + 1,
				     best_rmsd_array + myid*top_rmsd +j, chunk * sizeof (double));
                             memmove(best_triple_x_array[myid*top_rmsd + j + 1],
				     best_triple_x_array[myid*top_rmsd + j], chunk * 3 * sizeof (int));
                             memmove(best_triple_y_array[myid*top_rmsd + j + 1],
				     best_triple_y_array[myid*top_rmsd + j], chunk * 3 * sizeof (int));
                             memmove(best_quat_array[myid*top_rmsd + j + 1],
				     best_quat_array[myid*top_rmsd + j], chunk * 4 * sizeof (double));
                         }
                         best_rmsd_array[myid*top_rmsd + j] = best_rmsd_array[myid*top_rmsd +
									      stride *top_rmsd  + k];
                         memcpy(best_triple_x_array[myid*top_rmsd + j],
				best_triple_x_array[myid*top_rmsd + stride * top_rmsd + k], 3 * sizeof (int));
                         memcpy(best_triple_y_array[myid*top_rmsd + j],
				best_triple_y_array[myid*top_rmsd + stride * top_rmsd + k], 3 * sizeof (int));
                         memcpy(best_quat_array[myid*top_rmsd + j], best_quat_array[myid*top_rmsd + stride * top_rmsd + k], 4 * sizeof (double));
                         j++;
                         break; 
                     }
                     j++;
                 }
                 if (j == top_rmsd -1) break; // there is not any lower value
                 
             }
         }
     }
     
# endif     
     return 0;

 }



int find_best_triples_greedy(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
        double **best_quat) {

    int x_ctr, y_ctr, top_ctr;
    double **x = X_rep->full;
    int * x_type = X_rep->full_type;
    int NX = X_rep->N_full;
    double **y = Y_rep->full;
    int * y_type = Y_rep->full_type;
    int NY = Y_rep->N_full;
    int x_triple[3], y_triple[3];
    int chunk;
    double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
    double rmsd;
    double q_init[4] = {0.0};
    int done = 0;

    /***************************************/
    /* find reasonable triples of SSEs      */
    /* that correspond in type             */
    /*  and can be mapped onto each other  */
    /***************************************/

    for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
        best_rmsd[top_ctr] = BAD_RMSD + 1;
        best_triple_x[top_ctr][0] = -1;
    }




    for (x_ctr = 0; x_ctr < NX - 2 && !done; x_ctr++) {

        for (y_ctr = 0; y_ctr < NY - 2 && !done; y_ctr++) {

            if (y_type[y_ctr] != x_type[x_ctr]) continue;

            x_triple[0] = x_ctr;
            y_triple[0] = y_ctr;

            if (find_next_triple(x, y, x_type, y_type,
                    NX, NY, x_triple, y_triple)) {
                continue;
            }

            if (x_triple[1] < 0 || x_triple[2] < 0) continue;
            if (y_triple[1] < 0 || y_triple[2] < 0) continue;

            // do these three have  kind-of similar layout in space?
            // is handedness the same? 
            if (!same_hand_triple(X_rep, x_triple, Y_rep, y_triple, 3)) continue;
	    

            // are distances comparab;e? 
            if (distance_of_nearest_approach(X_rep, x_triple,
                    Y_rep, y_triple, 3, &rmsd)) continue;
            if (rmsd > cutoff_rmsd) continue;

            // find q_init that maps the two triples as well as possible
            if (opt_quat(x, NX, x_triple, y, NY, y_triple, 3, q_init, &rmsd)) continue;

            for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {

                if (rmsd <= best_rmsd[top_ctr]) {

                    chunk = no_top_rmsd - top_ctr - 1;

                    if (chunk) {
                        memmove(best_rmsd + top_ctr + 1, best_rmsd + top_ctr, chunk * sizeof (double));
                        memmove(best_quat[top_ctr + 1],
                                best_quat[top_ctr], chunk * 4 * sizeof (double));
                        memmove(best_triple_x[top_ctr + 1],
                                best_triple_x[top_ctr], chunk * 3 * sizeof (int));
                        memmove(best_triple_y[top_ctr + 1],
                                best_triple_y[top_ctr], chunk * 3 * sizeof (int));
                    }
                    best_rmsd[top_ctr] = rmsd;
                    memcpy(best_quat[top_ctr], q_init, 4 * sizeof (double));
                    memcpy(best_triple_x[top_ctr], x_triple, 3 * sizeof (int));
                    memcpy(best_triple_y[top_ctr], y_triple, 3 * sizeof (int));

                    break;
                }
            }

        }
    }
    return 0;
}


int find_submaps(int NX, int NY, Map * map, int * map_best){

    Map combined_map = {0};
    if (NX && NY) if ( initialize_map (&combined_map, NX, NY) ) exit (1);
    int i, j;
    int best_ctr = 0;
   
    while (  map_best[best_ctr] >  -1 ) {
	best_ctr ++;
    }

    
    if (best_ctr) {
	int nr_maps = (best_ctr<options.number_maps_cpl)?
	    best_ctr : options.number_maps_cpl;
	int best_i;
	int consistent;
	double z;
	double total_assigned_score, score, best_score = -100;
	double gap_score;

	for (i=0; i<nr_maps; i++) { /* look for the complement */
	    best_i =  map_best[i];
	    
	    /*intialize the (list of) submatch map(s) */
	    if ( !map[best_i].submatch_best) {
		/* for now look for a single map only */
		/* TODO - would it be worth any to look at more maps?*/ 
		int map_max = 1;
		map[best_i].submatch_best = emalloc (map_max*sizeof(int) );
		if (! map[best_i].submatch_best) return 1;
	    }
	    map[best_i].submatch_best[0]    = -1;
	    map[best_i].score_with_children =  0;
	    map[best_i].compl_z_score       =  0;
	    
	    for (j=0; j<best_ctr; j++) {

		if (i==j) continue;
		
		map_complementarity (map+best_i, map + map_best[j], &z);
			
		map_consistence (NX, NY, &combined_map, map+best_i, map + map_best[j],
				 &total_assigned_score, &gap_score);
		consistent = ((map+best_i)->assigned_score < total_assigned_score
			      && (map + map_best[j])->assigned_score < total_assigned_score);
		if ( consistent ) {
		    score = total_assigned_score;
		    if (  score > best_score ) {
			best_score = score;
			map[best_i].submatch_best[0] = map_best[j];
			map[best_i].score_with_children = total_assigned_score;
			map[best_i].compl_z_score = z;
		    }
		}
	    }

	}
    }
    
    free_map (&combined_map);
    return 0;
}
