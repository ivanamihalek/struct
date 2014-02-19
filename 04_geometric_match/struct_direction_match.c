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

# define  TOP_RMSD   200
# define  BAD_RMSD    10.0
# define  JACKFRUIT    8
# define  NUM_THREADS  8


int find_best_triples_exhaustive (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
				  double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
				  double **best_quat);
   
int find_best_triples_greedy(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
			     double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
			     double **best_quat);

int find_best_triples_exhaustive_parallel (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y, double **best_quat);
int find_best_triples_exhaustive_redux       (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
				  double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
					double **best_quat);

int find_best_triples_exhaustive_parallel_gpu (Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
				  double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
					double **best_quat);

int opt_quat (double ** x, int NX, int *set_of_directions_x,
	      double ** y, int NY, int *set_of_directions_y,
	      int set_size, double * q, double * rmsd){
    
    int opt_quat_lapack (double ** x, int NX, int *set_of_directions_x,
		  double ** y, int NY, int *set_of_directions_y,
		  int set_size, double * q, double * rmsd);
    int opt_quat_sine_lapack (double ** x, int NX, int *set_of_directions_x,
		  double ** y, int NY, int *set_of_directions_y,
		  int set_size, double * q, double * rmsd);
    
    return  opt_quat_sine_lapack (x, NX, set_of_directions_x,
			 y,  NY, set_of_directions_y,
			 set_size, q, rmsd);
}

/****************************************/
int direction_match (Representation* X_rep, Representation* Y_rep, List_of_maps *list){
	
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
    double **x_rotated    = NULL;
    double **tr_x_rotated = NULL;
    double **R;
    double z_scr = 0.0, *z_best;
    double avg, avg_sq, stdev;
    double alpha = options.alpha;
    double rmsd, *best_rmsd;
    double **best_quat;
    double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
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
    int stored_new;
    int * x2y, map_unstable;
    //time_t  time_now, time_start;
    Map * map     = list->map;
    int map_max   = list->no_maps_allocated;
    int *map_best = list->map_best;
   
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
    int qmap (double *x0, double *x1, double *y0, double *y1, double * quat);
    int store_sorted (Map * map, int NX, int NY, int *map_best, int map_max,
		      double * z_best, int best_ctr,
		      double z_scr, int  my_map_ctr, int *stored);
  
  

    smaller = (NX <= NY) ? NX : NY;
    no_top_rmsd = NX*NY/10; /* I'm not sure that this is the scale, but it works for now */

    /***********************/
    /* memory allocation   */
    /***********************/
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    if ( ! (x_rotated     = dmatrix (NX,3)) ) return 1;
    if ( ! (tr_x_rotated  = dmatrix (NX,3)) ) return 1;
    if ( ! (best_quat     = dmatrix (no_top_rmsd,4)) ) return 1;
    if ( ! (best_rmsd     = emalloc (no_top_rmsd*sizeof(double))) ) return 1;
    if ( ! (best_triple_x = intmatrix (no_top_rmsd,3)) ) return 1;
    if ( ! (best_triple_y = intmatrix (no_top_rmsd,3)) ) return 1;
    if ( ! (z_best = emalloc(NX*NY*sizeof(double) )) ) return 1;
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
    /* find reasonable triples of SSEs      */
    /* that correspond in type             */
    /*  and can be mapped onto each other  */
    /***************************************/


    //ProfilerStart("profile.out") ;
    /* the output will got to the file called profile.out */
    if (options.exhaustive) {
	/*
	 * Exhaustive search for all triplets
	 */
	if (options.gpu) {
# ifdef GPU
	    find_best_triples_exhaustive_parallel_gpu (X_rep, Y_rep, no_top_rmsd, best_rmsd, 
						       best_triple_x, best_triple_y, best_quat);
# else
	    fprintf (stderr, "from %s:%d: to use GPU pll, please recompile with -DGPU flag.\n",
		     __FILE__, __LINE__);
	    exit (1);
# endif
	  
	} else if (options.omp) {
# ifdef OMP
	    find_best_triples_exhaustive_parallel (X_rep, Y_rep, no_top_rmsd, best_rmsd, 
						   best_triple_x, best_triple_y, best_quat);
# else
	    fprintf (stderr, "from %s:%d: to use omp, please recompile with -DOMP flag.\n",
		     __FILE__, __LINE__);
	    exit (1);
# endif
	} else {
	    find_best_triples_exhaustive_redux (X_rep, Y_rep, no_top_rmsd, best_rmsd, 
						best_triple_x, best_triple_y, best_quat);
	}
    } else {
    
	/*
	 * Greedy search - old algorithm
	 */
  	find_best_triples_greedy (X_rep, Y_rep, no_top_rmsd, best_rmsd,  
				  best_triple_x, best_triple_y, best_quat);
    }
    
    //ProfilerStop();
    /*********************************************/
    /*   main loop                               */
    /*********************************************/
    map_ctr = 0;
    
    for (top_ctr=0; top_ctr<no_top_rmsd && done==0; top_ctr++) {
	if ( best_rmsd[top_ctr] > BAD_RMSD ) break;

	quat_to_R (best_quat[top_ctr], R);
	rotate (x_rotated, NX, R, x);

	F_current = F( y, y_type, NY, x_rotated, x_type, NX, alpha);

	/* find map which uses the 2 triples as anchors */
	no_anchors = 3;
	find_map (&penalty_params, X_rep, Y_rep, R, alpha, &F_effective, map + map_ctr,
		  best_triple_x[top_ctr], best_triple_y[top_ctr], no_anchors);

	/* does this map still map the two triples we started with? */
	x2y = (map + map_ctr) ->x2y;
	map_unstable  = 0;
	for (t=0; t<3; t++ ) {
	    if ( x2y[best_triple_x[top_ctr][t]] != best_triple_y[top_ctr][t] ) {
		map_unstable = 1;
	    }
	}
	if ( map_unstable) continue;
        /* do the mapped SSEs match in length? */
	if (options.use_length &&
	    (map+map_ctr)->avg_length_mismatch  > options.avg_length_mismatch_tol)  continue;
	
	
	/* dna here is not DNA but "distance of nearest approach" */
	cull_by_dna ( X_rep, best_triple_x[top_ctr], 
		      Y_rep, best_triple_y[top_ctr],  3,  map + map_ctr, cutoff_rmsd );
	

	/* monte that optimizes the aligned vectors only */
	for (i=0; i<NX; i++) {
	     x_type_fudg[i] = JACKFRUIT;
	}
	for (j=0; j<NY; j++) {
	     y_type_fudg[j] = JACKFRUIT*2;
	}

	no_anchors = 0;

	for (i=0; i<NX; i++) {
	     j = (map+map_ctr)->x2y[i];
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

	if (options.postprocess) {
	    z_scr = stdev ? (F_current - avg)/stdev : 0.0;
	} else {
	    z_scr =  0.0;
	}
	quat_to_R (q, R);
	/* store_sse_pair_score() is waste of time, but perhaps not critical */
	store_sse_pair_score (X_rep, Y_rep, R,  alpha, map + map_ctr);
	map_assigned_score (X_rep, map + map_ctr);


        /*   store the map that passed all the filters down to here*/
	map[map_ctr].F       = F_current;
	map[map_ctr].avg     = avg;
	map[map_ctr].avg_sq  = avg_sq;
	map[map_ctr].z_score = z_scr;
	memcpy ( map[map_ctr].q, q, 4*sizeof(double) );
		
	/* recalculate the assigned score*/

	/************************/
	/* store sorted         */
	/************************/
	/* find the place for the new z-score */
	store_sorted (map, NX, NY, map_best, map_max,
		      z_best, best_ctr, -map[map_ctr].assigned_score, map_ctr, &stored_new);

	if ( stored_new ) { /* we want to keep this map */
	    map_ctr ++;
	    best_ctr++;
	} /* otherwise this map space is reusable */
	    

	/* is this pretty much as good as it can get ? */
	if ( fabs (map[map_ctr].assigned_score - smaller)
	     < options.tol )  done = 1;

    }
    list->no_maps_used     = map_ctr;
    list->best_array_used  = best_ctr;

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
    free (z_best);
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
int store_sorted (Map * map, int NX, int NY, int *map_best,
		  int map_max,  double * z_best, int best_ctr,
		  double z_scr, int  map_ctr, int *stored) {

    int ctr, chunk;
    int correlated;
    double  z;
    
    *stored = 1; /* unless correlated and smaller in z,
		    we will store the map */
    
    /* is this map correlated with something
       that we have already ?*/
    correlated = 0;
    for  (ctr = 0; ctr<best_ctr && ! correlated && ctr < map_max; ctr ++ ) {
	map_complementarity ( map+map_ctr, map+ctr,  &z);
	correlated = (z < options.z_max_corr);
	/* TODO store calculated correlations, so I
	   don't have to redo it */
    }
    if ( ctr == map_max ) *stored=0;
    /* correlated: if this score is bigger, replace the old map */
    if ( correlated ) {
	if ( z_best[ctr] <  z_scr ) {
	      /* move */
	    chunk = best_ctr - ctr;
	    if ( chunk ) {
		memmove (z_best+ctr, z_best+ctr+1, chunk*sizeof(double)); 
		memmove (map_best+ctr, map_best+ctr+1, chunk*sizeof(int)); 
	    }
	    best_ctr --;
	} else {
	    /* we'll discard this map entirely */
	    *stored = 0;
	}
    }

    if ( *stored ) { /* if this map is to be stored */
	/* find the place in the list*/
	ctr = 0;
	while (ctr<best_ctr && ( z_scr >= z_best[ctr]) && ctr<map_max) ctr++;

	if ( ctr==map_max ) {/* we alredy have enough maps
			       which are much better */
	    *stored = 0;
	} else {
	    /* move */
	    chunk = best_ctr - ctr;
	    if ( chunk ) {
		memmove (z_best+ctr+1, z_best+ctr, chunk*sizeof(double)); 
		memmove (map_best+ctr+1, map_best+ctr, chunk*sizeof(int)); 
	    }
	    /* store z_score*/
	    z_best[ctr]   = z_scr;
	    map_best[ctr] = map_ctr;
	}
    }   

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
    /* numerically negative? algebraically negative? toplogically negative? epistemologically negative? */
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
int opt_quat_old ( double ** x, int NX, int *set_of_directions_x,
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
    int chunk;
    int NX = X_rep->N_full;
    int NY = Y_rep->N_full;
    int * x_type = X_rep->full_type;
    int * y_type = Y_rep->full_type;
    TripleID *x_triple, *y_triple;
 
    double **x = X_rep->full;    
    double **y = Y_rep->full;
    
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

    if ( ! (x_triple = emalloc (MAX_TRIPS*sizeof(TripleID) ) ) )  return 1;
    if ( ! (y_triple = emalloc (MAX_TRIPS*sizeof(TripleID) ) ) )  return 1;
    
    no_trip_pairs_to_compare = 0;
    x_enumeration_done       = 0;
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
	    /* if we are doing out-of-order search we should try all permutations of these triplets
	     storage space issues? */
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

	    if ( hand(X_rep, x_triple[xtrip_ct].member ) ) x_triple[xtrip_ct].fingerprint |= HAND;
		
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

    exit (1);
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
    int chunk;
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
