# include "struct.h"

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
       int he angle makeing them pointing toward each other or
       away from each other changes the direction of cross prod;
       rather, use one vector, and distance between the cm's as the other */
    for (i=0; i<3; i++ ) {
	cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
    }
    normalized_cross (x[ray_a], x[ray_b], cross, &aux); 
    /* note I am making another cm vector here */
    for (i=0; i<3; i++ ) {
	avg_cm[i] = (x_cm[ray_b][i] + x_cm[ray_a][i])/2;
	cm_vector[i] = x_cm[ray_c][i] - avg_cm[i];
    }
    unnorm_dot (cm_vector, cross, &d);

    if ( d > 0 ) {
	return 0;
    } else {
	return 1;
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
	    avg_cm[i]    = (x_cm[ray_b][i] + x_cm[ray_a][i])/2;
	    cm_vector[i] =  x_cm[ray_c][i] - avg_cm[i];
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
	if (normalized_cross (y[ray_a], cm_vector, cross, &aux)) {
	    /* if  I can calculate the cross product for one triplet
	       but not for the other, I declare it a different hand (this will
	       be dropped from further consideration) */
	    return 0;
	}
	/* note I am making another cm vector here */
	for (i=0; i<3; i++ ) {
	    avg_cm[i]    = (y_cm[ray_b][i] + y_cm[ray_a][i])/2;
	    cm_vector[i] =  y_cm[ray_c][i] - avg_cm[i];
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

	/* if the vectors linking cms are parallel, the hand is the same (my defintion) */
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
	    /* distance of nearest approach of ray a
	       to the cm of b, in the set of directions x */
	    for (i=0; i<3; i++ ) {
		cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
	    }
	    /* the norm of the cross product is ||cm_vector||*sin(alpha),
	       bc the norm of  the rep vectors is 1 */
	    normalized_cross (cm_vector, x[ray_a], cross, &distance_x);

	    ray_a = set_of_directions_y[a];
	    ray_b = set_of_directions_y[b];
 	    /* distance of nearest approach of ray a
	       to the cm of b, in the set of directions y */
	    for (i=0; i<3; i++ ) {
		cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
	    }
	    normalized_cross (cm_vector, y[ray_a], cross, &distance_y);


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
    
    //printf ("***** rmsd:  %8.3lf  \n\n", rmsd);

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
    int i, j, ray_a, ray_b, a;
    int in_triple, coord_ctr, norm;
    
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
	/* the rmsd for the dna of this element to the triple */
	for (a=0; a<set_size; a++) {

	    /**************************************/
	    ray_a = set_of_directions_x[a];
	    ray_b = i;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions x */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = x_cm[ray_b][coord_ctr] - x_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, x[ray_a], cross, &distance_x);

	    ray_a = set_of_directions_y[a];
	    ray_b = j;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions y */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = y_cm[ray_b][coord_ctr] - y_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, y[ray_a], cross, &distance_y);

	    aux   = distance_x-distance_y;
	    rmsd += aux*aux;
	    norm ++;

	    /**************************************/
	    ray_b = set_of_directions_x[a];
	    ray_a = i;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions x */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = x_cm[ray_b][coord_ctr] - x_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, x[ray_a], cross, &distance_x);

	    ray_b = set_of_directions_y[a];
	    ray_a = j;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions y */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = y_cm[ray_b][coord_ctr] - y_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, y[ray_a], cross, &distance_y);

	    aux   = distance_x-distance_y;
	    rmsd += aux*aux;
	    norm ++;

	}

	rmsd /= norm;
	rmsd = sqrt(rmsd);

	/* if rmsd is bigger than the cutoff, lose this from the mapping */
	if ( rmsd > cutoff_rmsd) {
	    map->x2y[i] = -1;
	    map->y2x[j] = -1;
	}
    }
    
     
    return 0;
}

/********************************************************/
int init_triples_array(Triples_array * triples_array, int size){
     if ( ! (triples_array->hhh_array  = intmatrix(size,3)) ) return 0;
     if ( ! (triples_array->hhs_array  = intmatrix(size,3)) ) return 0;
     if ( ! (triples_array->hsh_array  = intmatrix(size,3)) ) return 0;
     if ( ! (triples_array->hss_array  = intmatrix(size,3)) ) return 0;
     if ( ! (triples_array->shh_array  = intmatrix(size,3)) ) return 0;
     if ( ! (triples_array->ssh_array  = intmatrix(size,3)) ) return 0;
     if ( ! (triples_array->shs_array  = intmatrix(size,3)) ) return 0;
     if ( ! (triples_array->sss_array  = intmatrix(size,3)) ) return 0;
     
     return 1;
}
/********************************************************/
int free_triples_array(Triples_array * triples_array) {
    free_imatrix(triples_array->hhh_array); 
    free_imatrix(triples_array->hhs_array);
    free_imatrix(triples_array->hsh_array);
    free_imatrix(triples_array->hss_array);
    free_imatrix(triples_array->shh_array);
    free_imatrix(triples_array->ssh_array); 
    free_imatrix(triples_array->shs_array); 
    free_imatrix(triples_array->sss_array);
    return 0;
}

    
/********************************************************/
/********************************************************/
/********************************************************/
int sortTriplets(int ** best_triple_x_array, int ** best_triple_y_array,
		 double * best_rmsd_array, double ** best_quat_array, int top_rmsd) {

    int j, k, chunk;
    
    j = 0;
    for (k = 0; k < top_rmsd; ++k) {
	while (j < top_rmsd) {
	    if (best_rmsd_array[k] < best_rmsd_array[j]) {

		chunk = top_rmsd - j - 1;

		if (chunk) {
		    memmove(best_rmsd_array +  j + 1,
			    best_rmsd_array +  j, chunk * sizeof (double));
		    memmove(best_triple_x_array[ j + 1],
			    best_triple_x_array[ j], chunk * 3 * sizeof (int));
		    memmove(best_triple_y_array[ j + 1],
			    best_triple_y_array[ j], chunk * 3 * sizeof (int));
		}
		best_rmsd_array[ j] = best_rmsd_array[k];
		memcpy(best_triple_x_array[ j],
		       best_triple_x_array[ k], 3 * sizeof (int));
		memcpy(best_triple_y_array[ j],
		       best_triple_y_array[ k], 3 * sizeof (int));
		j++;
		break;
	    }
	    j++;
	}
	if (j == top_rmsd - 1) break; // there is not any lower value

    }
    return 0;

}

// function that inserts triples to an array of triples using types of SSE in triples as a key

int insert_triple_to_array(int * triple_types, int triple1, int triple2, int triple3, Triples_array *triple_array) {
    int key;

    key = triple_types[triple1] * 100 + triple_types[triple2] * 10 + triple_types[triple3];

    // s eq 4, h eq 2    
    switch (key) {
        case 222:
            triple_array->hhh_array[triple_array->hhh_cnt][0] = triple1;
            triple_array->hhh_array[triple_array->hhh_cnt][1] = triple2;
            triple_array->hhh_array[triple_array->hhh_cnt][2] = triple3;
            triple_array->hhh_cnt++;
            break;
        case 224:
            triple_array->hhs_array[triple_array->hhs_cnt][0] = triple1;
            triple_array->hhs_array[triple_array->hhs_cnt][1] = triple2;
            triple_array->hhs_array[triple_array->hhs_cnt][2] = triple3;
            triple_array->hhs_cnt++;
            break;
        case 242:
            triple_array->hsh_array[triple_array->hsh_cnt][0] = triple1;
            triple_array->hsh_array[triple_array->hsh_cnt][1] = triple2;
            triple_array->hsh_array[triple_array->hsh_cnt][2] = triple3;
            triple_array->hsh_cnt++;
            break;
        case 244:
            triple_array->hss_array[triple_array->hss_cnt][0] = triple1;
            triple_array->hss_array[triple_array->hss_cnt][1] = triple2;
            triple_array->hss_array[triple_array->hss_cnt][2] = triple3;
            triple_array->hss_cnt++;
            break;
        case 422:
            triple_array->shh_array[triple_array->shh_cnt][0] = triple1;
            triple_array->shh_array[triple_array->shh_cnt][1] = triple2;
            triple_array->shh_array[triple_array->shh_cnt][2] = triple3;
            triple_array->shh_cnt++;
            break;
        case 442:
            triple_array->ssh_array[triple_array->ssh_cnt][0] = triple1;
            triple_array->ssh_array[triple_array->ssh_cnt][1] = triple2;
            triple_array->ssh_array[triple_array->ssh_cnt][2] = triple3;
            triple_array->ssh_cnt++;
            break;
        case 424:
            triple_array->shs_array[triple_array->shs_cnt][0] = triple1;
            triple_array->shs_array[triple_array->shs_cnt][1] = triple2;
            triple_array->shs_array[triple_array->shs_cnt][2] = triple3;
            triple_array->shs_cnt++;
            break;
        case 444:
            triple_array->sss_array[triple_array->sss_cnt][0] = triple1;
            triple_array->sss_array[triple_array->sss_cnt][1] = triple2;
            triple_array->sss_array[triple_array->sss_cnt][2] = triple3;
            triple_array->sss_cnt++;
            break;
        default:
            return 0;
    }

    return 1;

}

int insert_triple_to_heap(Representation* X_rep, Representation* Y_rep,
			  int ** x_triple_array, int ** y_triple_array, int x_triple_cnt, int y_triple_cnt, PriorityQueue * heap) {
    

    Triple * triple_array = malloc(y_triple_cnt * sizeof(Triple));
    int m;
    for (m = 0; m < y_triple_cnt; ++m) triple_array[m].rmsd = 11;
    
    int i, j;
    int cnt_x = x_triple_cnt;
    int cnt_y = y_triple_cnt;
    
    int ** x_triple_array_local = intmatrix(cnt_x, 3);
    int ** y_triple_array_local = intmatrix(cnt_y, 3);
            
    memcpy(*x_triple_array_local, *x_triple_array, sizeof (int) * 3 * cnt_x);
    memcpy(*y_triple_array_local, *y_triple_array, sizeof (int) * 3 * cnt_y);
    
    int x_triple[3];
    int y_triple[3];
    double rmsd;
    double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
    double q_init[4] = {0.0}; // no change
    
    for (i = 0; i < cnt_x; ++i) {
        for (j = 0; j < cnt_y; ++j) {
            x_triple[0] = x_triple_array_local[i][0];
            x_triple[1] = x_triple_array_local[i][1];
            x_triple[2] = x_triple_array_local[i][2];
    
            
            y_triple[0] = y_triple_array_local[j][0];
            y_triple[1] = y_triple_array_local[j][1];
            y_triple[2] = y_triple_array_local[j][2];

	    if (!same_hand_triple(X_rep, x_triple, Y_rep, y_triple, 3)) continue;
	    if (distance_of_nearest_approach(X_rep, x_triple, Y_rep, y_triple, 3, &rmsd)) continue;
	    if (rmsd > cutoff_rmsd) continue;

            triple_array[j].rmsd = rmsd;
            memcpy(&triple_array[j].quat, &q_init, 4 * sizeof (double));
            memcpy(&triple_array[j].triple_x, &x_triple, 3 * sizeof (double));
            memcpy(&triple_array[j].triple_y, &y_triple, 3 * sizeof (double));
        }
    } 

    free(triple_array);
    
    return 0;
} 

