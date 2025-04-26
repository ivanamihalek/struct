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
# include "struct_triplets.h"
# include "sys/time.h"
# ifdef OMP
#   include "omp.h"
# endif

//# include "gperftools/profiler.h"

# define  TOP_RMSD   200
# define  BAD_RMSD    10.0
# define  JACKFRUIT    8
# define  NUM_THREADS  8

int find_best_triples_exhaustive(Representation* X_rep, Representation* Y_rep,
		int no_top_rmsd, double * best_rmsd, int ** best_triple_x,
		int ** best_triple_y, double **best_quat);

int find_best_triples_greedy(Representation* X_rep, Representation* Y_rep,
		int no_top_rmsd, double * best_rmsd, int ** best_triple_x,
		int ** best_triple_y, double **best_quat);

int find_best_triples_exhaustive_parallel(Representation* X_rep,
		Representation* Y_rep, int no_top_rmsd, double * best_rmsd,
		int ** best_triple_x, int ** best_triple_y, double **best_quat);
int find_best_triples_exhaustive_redux(Representation* X_rep,
		Representation* Y_rep, int no_top_rmsd, double * best_rmsd,
		int ** best_triple_x, int ** best_triple_y, double **best_quat);

int find_best_triples_exhaustive_parallel_gpu(Representation* X_rep,
		Representation* Y_rep, int no_top_rmsd, double * best_rmsd,
		int ** best_triple_x, int ** best_triple_y, double **best_quat);

int opt_quat(double ** x, int NX, int *set_of_directions_x, double ** y, int NY,
		int *set_of_directions_y, int set_size, double * q, double * rmsd);

/****************************************/
int direction_match(Representation* X_rep, Representation* Y_rep,
		List_of_maps *list) {

	Penalty_parametrization penalty_params; /* for SW */
	double **x = X_rep->full;
	int * x_type = X_rep->full_type;
	int NX = X_rep->N_full;
	double **y = Y_rep->full;
	int * y_type = Y_rep->full_type;
	int NY = Y_rep->N_full;

	double F_effective = 0.0;
	double F_current;
	double q[4] = { 0.0 };
	double **x_rotated = NULL;
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
	int map_ctr = 0;
	int stored_new;
	int * x2y, map_unstable;
	//time_t  time_now, time_start;
	Map * map = list->map;
	int map_max = list->no_maps_allocated;
	int *map_best = list->map_best;

	int cull_by_dna(Representation * X_rep, int *set_of_directions_x,
			Representation * Y_rep, int *set_of_directions_y, int set_size,
			Map *map, double cutoff_rmsd);
	int distance_of_nearest_approach(Representation * X_rep,
			int *set_of_directions_x, Representation * Y_rep,
			int *set_of_directions_y, int set_size, double * rmsd_ptr);
	int same_hand_triple(Representation * X_rep, int *set_of_directions_x,
			Representation * Y_rep, int *set_of_directions_y, int set_size);

	int find_map(Penalty_parametrization * params, Representation *X_rep,
			Representation *Y_rep, double ** R, double alpha,
			double * F_effective, Map *map, int *anchor_x, int * anchor_y,
			int anchor_size);
	int find_next_triple(double **X, double **Y, int *x_type, int *y_type,
			int NX, int NY, int *x_triple, int *y_triple);
	int gradient_descent(int first_call, double alpha, double **x, int *x_type,
			int NX, double **y, int *y_type, int NY, double *q_best,
			double *F_best_ptr);
	int monte_carlo(double alpha, double **x, int * x_type, int NX, double **y,
			int * y_type, int NY, double *q_best, double *F_best_ptr);
	int qmap(double *x0, double *x1, double *y0, double *y1, double * quat);
	int store_sorted(Map * map, int NX, int NY, int *map_best, int map_max,
			double * z_best, int best_ctr, double z_scr, int my_map_ctr,
			int *stored);

	smaller = (NX <= NY) ? NX : NY;
	no_top_rmsd = NX * NY / 10; /* I'm not sure that this is the scale, but it works for now */
	if (options.current_algorithm != SEQUENTIAL) {
		no_top_rmsd *= 10;
	}

	/***********************/
	/* memory allocation   */
	/***********************/
	if (!(R = dmatrix(3, 3)))
		return 1; /* compiler is bugging me otherwise */
	if (!(x_rotated = dmatrix(NX, 3)))
		return 1;
	if (!(tr_x_rotated = dmatrix(NX, 3)))
		return 1;
	if (!(best_quat = dmatrix(no_top_rmsd, 4)))
		return 1;
	if (!(best_rmsd = emalloc(no_top_rmsd * sizeof(double))))
		return 1;
	if (!(best_triple_x = intmatrix(no_top_rmsd, 3)))
		return 1;
	if (!(best_triple_y = intmatrix(no_top_rmsd, 3)))
		return 1;
	if (!(z_best = emalloc(NX * NY * sizeof(double))))
		return 1;
	if (!(x_type_fudg = emalloc(NX * sizeof(int))))
		return 1;
	if (!(y_type_fudg = emalloc(NY * sizeof(int))))
		return 1;
	if (!(anchor_x = emalloc(NX * sizeof(int))))
		return 1;
	if (!(anchor_y = emalloc(NY * sizeof(int))))
		return 1;
	/***********************/

	/***********************/
	/* expected quantities */
	/***********************/
	avg = avg_sq = stdev = 0.0;
	if (0) {
	    if (!options.path[0]) {
		 fprintf(stderr,  "path to integral table must be given "
				  "in the cmd file (kwd \"path\") to calculate z-score for F.\n");
		 exit(1);
	     }
	     if (F_moments(x, x_type, NX, y, y_type, NY, alpha, &avg, &avg_sq, &stdev))
			return 1;
	}
	/***********************/

	/***********************/
	/* initialization      */
	/***********************/
	best_ctr = 0;
	penalty_params.gap_opening = options.gap_open;
	penalty_params.gap_extension = options.gap_extend;
	penalty_params.endgap = options.endgap;
	penalty_params.endgap_special_treatment = options.use_endgap;
	penalty_params.custom_gap_penalty_x = NULL;
	penalty_params.custom_gap_penalty_y = NULL;
	if (!(penalty_params.custom_gap_penalty_x = emalloc(NX * sizeof(double))))
		return 1;
	if (!(penalty_params.custom_gap_penalty_y = emalloc(NY * sizeof(double))))
		return 1;
	/***********************/

	/***************************************/
	/* find reasonable triples of SSEs      */
	/* that correspond in type             */
	/*  and can be mapped onto each other  */
	/***************************************/

	//ProfilerStart("profile.out") ;
	/* the output will got to the file called profile.out */
	find_best_triples_exhaustive_redux(X_rep, Y_rep, no_top_rmsd, best_rmsd,
			best_triple_x, best_triple_y, best_quat);

	//ProfilerStop();
	/*********************************************/
	/*   main loop                               */
	/*********************************************/
	map_ctr = 0;

	for (top_ctr = 0; top_ctr < no_top_rmsd && done == 0; top_ctr++) {

		if (best_rmsd[top_ctr] > BAD_RMSD) break;

		quat_to_R(best_quat[top_ctr], R);
		rotate(x_rotated, NX, R, x);

		F_current = F(y, y_type, NY, x_rotated, x_type, NX, alpha);

		/* find map which uses the 2 triples as anchors */
		no_anchors = 3;
		find_map(&penalty_params, X_rep, Y_rep, R, alpha, &F_effective,
				map + map_ctr, best_triple_x[top_ctr], best_triple_y[top_ctr],
				no_anchors);

		/* does this map still map the two triples we started with? */
		x2y = (map + map_ctr)->x2y;
		map_unstable = 0;
		for (t = 0; t < 3; t++) {
			if (x2y[best_triple_x[top_ctr][t]] != best_triple_y[top_ctr][t]) {
				map_unstable = 1;
				break;
			}
		}
		if (map_unstable) continue;

		/* do the mapped SSEs match in length? */
		if (options.use_length && (map + map_ctr)->avg_length_mismatch > options.avg_length_mismatch_tol)
			continue;

		/* dna here is not DNA but "distance of nearest approach" */
		cull_by_dna(X_rep, best_triple_x[top_ctr], Y_rep,
				best_triple_y[top_ctr], 3, map + map_ctr, cutoff_rmsd);

		/* monte that optimizes the aligned vectors only */
		for (i = 0; i < NX; i++) {
			x_type_fudg[i] = JACKFRUIT;
		}
		for (j = 0; j < NY; j++) {
			y_type_fudg[j] = JACKFRUIT * 2;
		}

		no_anchors = 0;

		for (i = 0; i < NX; i++) {
			j = (map + map_ctr)->x2y[i];
			if (j < 0)
				continue;
			x_type_fudg[i] = x_type[i];
			y_type_fudg[j] = y_type[j];
			anchor_x[no_anchors] = i;
			anchor_y[no_anchors] = j;
			no_anchors++;
		}

		if (opt_quat(x, NX, anchor_x, y, NY, anchor_y, no_anchors, q, &rmsd))
			continue;

		retval = monte_carlo(alpha, x, x_type_fudg, NX, y, y_type_fudg, NY, q,
				&F_current);
		if (retval)
			return retval;

		if (options.postprocess) {
			z_scr = stdev ? (F_current - avg) / stdev : 0.0;
		} else {
			z_scr = 0.0;
		}
		quat_to_R(q, R);
		/* store_sse_pair_score() is waste of time, but perhaps not critical */
		store_sse_pair_score(X_rep, Y_rep, R, alpha, map + map_ctr);
		map_assigned_score(X_rep, map + map_ctr);

		/*   store the map that passed all the filters down to here*/
		map[map_ctr].F = F_current;
		map[map_ctr].avg = avg;
		map[map_ctr].avg_sq = avg_sq;
		map[map_ctr].z_score = z_scr;
		memcpy(map[map_ctr].q, q, 4 * sizeof(double));

		/* recalculate the assigned score*/

		/************************/
		/* store sorted         */
		/************************/
		/* find the place for the new z-score */
		store_sorted(map, NX, NY, map_best, map_max, z_best, best_ctr,
				-map[map_ctr].assigned_score, map_ctr, &stored_new);

		if (stored_new) { /* we want to keep this map */
			map_ctr++;
			best_ctr++;
		} /* otherwise this map space is reusable */

		/* is this pretty much as good as it can get ?
		 if ( fabs (map[map_ctr].assigned_score - smaller)
		 < options.tol )  done = 1;
		 */
	}
	list->no_maps_used = map_ctr;
	list->best_array_used = best_ctr;

	/**********************/
	/* garbage collection */
	gradient_descent(1, 0.0, NULL, NULL, 0,
	NULL, NULL, 0, NULL, NULL);
	free_dmatrix(R);
	free_dmatrix(x_rotated);
	free_dmatrix(tr_x_rotated);
	free_dmatrix(best_quat);
	free_imatrix(best_triple_x);
	free_imatrix(best_triple_y);
	free(z_best);
	free(best_rmsd);
	free(x_type_fudg);
	free(y_type_fudg);
	free(anchor_x);
	free(anchor_y);

	if (penalty_params.custom_gap_penalty_x)
		free(penalty_params.custom_gap_penalty_x);
	if (penalty_params.custom_gap_penalty_y)
		free(penalty_params.custom_gap_penalty_y);
	/*********************/

	return 0;
}

/***************************************/
/***************************************/
int store_sorted(Map * map, int NX, int NY, int *map_best, int map_max,
		double * z_best, int best_ctr, double z_scr, int map_ctr, int *stored) {

	int ctr, chunk;
	int correlated;
	double z;

	*stored = 1; /* unless correlated and smaller in z,
	 we will store the map */

	/* is this map correlated with something
	 that we have already ?*/
	correlated = 0;
	for (ctr = 0; ctr < best_ctr && !correlated && ctr < map_max; ctr++) {
		map_complementarity(map + map_ctr, map + ctr, &z);
		correlated = (z < options.z_max_corr);
		/* TODO store calculated correlations, so I
		 don't have to redo it */
	}
	if (ctr == map_max)
		*stored = 0;
	/* correlated: if this score is bigger, replace the old map */
	if (correlated) {
		if (z_best[ctr] < z_scr) {
			/* move */
			chunk = best_ctr - ctr;
			if (chunk) {
				memmove(z_best + ctr, z_best + ctr + 1, chunk * sizeof(double));
				memmove(map_best + ctr, map_best + ctr + 1,
						chunk * sizeof(int));
			}
			best_ctr--;
		} else {
			/* we'll discard this map entirely */
			*stored = 0;
		}
	}

	if (*stored) { /* if this map is to be stored */
		/* find the place in the list*/
		ctr = 0;
		while (ctr < best_ctr && (z_scr >= z_best[ctr]) && ctr < map_max)
			ctr++;

		if (ctr == map_max) {/* we alredy have enough maps
		 which are much better */
			*stored = 0;
		} else {
			/* move */
			chunk = best_ctr - ctr;
			if (chunk) {
				memmove(z_best + ctr + 1, z_best + ctr, chunk * sizeof(double));
				memmove(map_best + ctr + 1, map_best + ctr,
						chunk * sizeof(int));
			}
			/* store z_score*/
			z_best[ctr] = z_scr;
			map_best[ctr] = map_ctr;
		}
	}

	return 0;
}

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

int find_best_triples_exhaustive_redux(Representation* X_rep,
		Representation* Y_rep, int no_top_rmsd, double * best_rmsd,
		int ** best_triple_x, int ** best_triple_y, double **best_quat) {

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
	double q_init[4] = { 0.0 };
	double threshold_dist = options.threshold_distance;

	double ** cmx = X_rep->cm;
	double ** cmy = Y_rep->cm;

	/* for out of order search */
	int number_of_permutations = 1;
	int permutation[6][3] = { { 0 } };
	int p;

	if (options.verbose)
		printf("exhaustive search \n");

	/***************************************/
	/* find reasonable triples of SSEs      */
	/* that correspond in type             */
	/*  and can be mapped onto each other  */
	/***************************************/
	for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
		best_rmsd[top_ctr] = BAD_RMSD + 1;
		best_triple_x[top_ctr][0] = -1;
	}

	/* max trips is currently set to 5000 in struct_triplets.h  */
	/* te hope is that there won't actually be that many triplets passing the
	 distance filter */
	if (!(x_triple = emalloc(MAX_TRIPS * sizeof(TripleID))))
		return 1;
	if (!(y_triple = emalloc(MAX_TRIPS * sizeof(TripleID))))
		return 1;

	no_trip_pairs_to_compare = 0;
	x_enumeration_done = 0;
	i_x = 0;
	j_x = 1;
	k_x = 1;

	panic_ctr = 0;
	while (!x_enumeration_done) {

		x_list_full = 0;
		xtrip_ct = 0;

		/***********************************************************************/
		while (!x_list_full && !x_enumeration_done) {

			k_x++;
			if (k_x == NX) {

				if (j_x < NX - 2) {
					j_x++;
					k_x = j_x + 1;

				} else {

					if (i_x < NX - 3) {
						i_x++;
						j_x = i_x + 1;
						k_x = j_x + 1;

					} else {
						x_enumeration_done = 1;
						break; /* from filling the x list */
					}
				}
			}

			if (two_point_distance(cmx[i_x], cmx[j_x]) > threshold_dist)
				continue;
			if (two_point_distance(cmx[i_x], cmx[k_x]) > threshold_dist)
				continue;
			if (two_point_distance(cmx[j_x], cmx[k_x]) > threshold_dist)
				continue;

			x_triple[xtrip_ct].fingerprint = 0;
			if (x_type[i_x] == HELIX)
				x_triple[xtrip_ct].fingerprint |= TYPE1;
			if (x_type[j_x] == HELIX)
				x_triple[xtrip_ct].fingerprint |= TYPE2;
			if (x_type[k_x] == HELIX)
				x_triple[xtrip_ct].fingerprint |= TYPE3;

			x_triple[xtrip_ct].member[0] = i_x;
			x_triple[xtrip_ct].member[1] = j_x;
			x_triple[xtrip_ct].member[2] = k_x;

			if (hand(X_rep, x_triple[xtrip_ct].member))
				x_triple[xtrip_ct].fingerprint |= HAND;

			xtrip_ct++;
			x_list_full = (xtrip_ct == MAX_TRIPS);

		} /* end filling the x trips list */

		no_xtrips = xtrip_ct;

		y_enumeration_done = 0;
		i_y = 0;
		j_y = 1;
		k_y = 1;
		/***********************************************************************/
		while (!y_enumeration_done) {

			y_list_full = 0;
			ytrip_ct = 0;

			while (!y_list_full && !y_enumeration_done) {

				k_y++;
				if (k_y == NY) {

					if (j_y < NY - 2) {
						j_y++;
						k_y = j_y + 1;

					} else {

						if (i_y < NY - 3) {
							i_y++;
							j_y = i_y + 1;
							k_y = j_y + 1;

						} else {
							y_enumeration_done = 1;
							break; /* from filling the y list */
						}
					}
				}

				if (two_point_distance(cmy[i_y], cmy[j_y]) > threshold_dist)
					continue;
				if (two_point_distance(cmy[i_y], cmy[k_y]) > threshold_dist)
					continue;
				if (two_point_distance(cmy[j_y], cmy[k_y]) > threshold_dist)
					continue;

				/* if we are doing out-of-order search we should try all
				 permutations of these triplets  - storage space issues? */

				if (options.current_algorithm == SEQUENTIAL) {
					number_of_permutations = 1;
					permutation[0][0] = i_y;
					permutation[0][1] = j_y;
					permutation[0][2] = k_y;
				} else {
					number_of_permutations = 6;
					permutation[0][0] = i_y;
					permutation[0][1] = j_y;
					permutation[0][2] = k_y;
					permutation[1][0] = i_y;
					permutation[1][1] = k_y;
					permutation[1][2] = j_y;
					permutation[2][0] = j_y;
					permutation[2][1] = k_y;
					permutation[2][2] = i_y;
					permutation[3][0] = j_y;
					permutation[3][1] = i_y;
					permutation[3][2] = k_y;
					permutation[4][0] = k_y;
					permutation[4][1] = i_y;
					permutation[4][2] = j_y;
					permutation[5][0] = k_y;
					permutation[5][1] = j_y;
					permutation[5][2] = i_y;
				}

				for (p = 0; p < number_of_permutations; p++) {

					int i_p = permutation[p][0];
					int j_p = permutation[p][1];
					int k_p = permutation[p][2];

					y_triple[ytrip_ct].fingerprint = 0;
					if (y_type[i_p] == HELIX)
						y_triple[ytrip_ct].fingerprint |= TYPE1;
					if (y_type[j_p] == HELIX)
						y_triple[ytrip_ct].fingerprint |= TYPE2;
					if (y_type[k_p] == HELIX)
						y_triple[ytrip_ct].fingerprint |= TYPE3;

					y_triple[ytrip_ct].member[0] = i_p;
					y_triple[ytrip_ct].member[1] = j_p;
					y_triple[ytrip_ct].member[2] = k_p;

					if (hand(Y_rep, y_triple[ytrip_ct].member))
						y_triple[ytrip_ct].fingerprint |= HAND;

					ytrip_ct++;
					y_list_full = (ytrip_ct == MAX_TRIPS);
					if (y_list_full)
						break;
				}

			} /* end filling the y trips list */

			no_ytrips = ytrip_ct;
			panic_ctr++;
			if (panic_ctr == 1000) {
				fprintf(stderr, "Exit by panic ctr in %s:%d.\n", __FILE__,
						__LINE__);
				exit(1);
			}

			/*****************************************************/
			for (xtrip_ct = 0; xtrip_ct < no_xtrips; xtrip_ct++) {
				for (ytrip_ct = 0; ytrip_ct < no_ytrips; ytrip_ct++) {

					no_trip_pairs_to_compare++;
					/* filters :*/
					if (x_triple[xtrip_ct].fingerprint
							^ y_triple[ytrip_ct].fingerprint)
						continue;

					if (distance_of_nearest_approach(X_rep,
							x_triple[xtrip_ct].member, Y_rep,
							y_triple[ytrip_ct].member, 3, &rmsd))
						continue;
					if (rmsd > cutoff_rmsd)
						continue;

					if (opt_quat(x, NX, x_triple[xtrip_ct].member, y, NY,
							y_triple[ytrip_ct].member, 3, q_init, &rmsd))
						continue;

					/* store the survivors */
					for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {

						if (rmsd <= best_rmsd[top_ctr]) {

							chunk = no_top_rmsd - top_ctr - 1;

							if (chunk) {
								memmove(best_rmsd + top_ctr + 1,
										best_rmsd + top_ctr,
										chunk * sizeof(double));
								memmove(best_quat[top_ctr + 1],
										best_quat[top_ctr],
										chunk * 4 * sizeof(double));
								memmove(best_triple_x[top_ctr + 1],
										best_triple_x[top_ctr],
										chunk * 3 * sizeof(int));
								memmove(best_triple_y[top_ctr + 1],
										best_triple_y[top_ctr],
										chunk * 3 * sizeof(int));
							}
							best_rmsd[top_ctr] = rmsd;
							memcpy(best_quat[top_ctr], q_init,
									4 * sizeof(double));
							memcpy(best_triple_x[top_ctr],
									x_triple[xtrip_ct].member, 3 * sizeof(int));
							memcpy(best_triple_y[top_ctr],
									y_triple[ytrip_ct].member, 3 * sizeof(int));

							break;

						}
					}

				}
			}
			/*****************************************************/

		} /* y enumeration loop */
	} /* x enumeration loop  */

# if 0
	printf ("no trips in x: %d\n", no_xtrips);
	printf ("no trips in y: %d\n", no_ytrips);
	printf ("no trips to compare : %d\n", no_trip_pairs_to_compare);

	for (top_ctr=0; top_ctr<no_top_rmsd; top_ctr++) {

		if ( best_rmsd[top_ctr] > BAD_RMSD ) break;
		printf (" %8.3lf      %3d  %3d  %3d     %3d  %3d  %3d \n", best_rmsd[top_ctr],
				best_triple_x[top_ctr][0], best_triple_x[top_ctr][1], best_triple_x[top_ctr][2],
				best_triple_y[top_ctr][0], best_triple_y[top_ctr][1], best_triple_y[top_ctr][2]);
	}

	// exit (1);
# endif

	free(x_triple);
	free(y_triple);

	return 0;

}

int find_submaps(int NX, int NY, Map * map, int * map_best) {

	Map combined_map = { 0 };
	if (NX && NY)
		if (initialize_map(&combined_map, NX, NY))
			exit(1);
	int i, j;
	int best_ctr = 0;

	while (map_best[best_ctr] > -1) {
		best_ctr++;
	}

	if (best_ctr) {
		int nr_maps =
				(best_ctr < options.number_maps_cpl) ?
						best_ctr : options.number_maps_cpl;
		int best_i;
		int consistent;
		double z;
		double total_assigned_score, score, best_score = -100;
		double gap_score;

		for (i = 0; i < nr_maps; i++) { /* look for the complement */
			best_i = map_best[i];

			/*intialize the (list of) submatch map(s) */
			if (!map[best_i].submatch_best) {
				/* for now look for a single map only */
				/* TODO - would it be worth any to look at more maps?*/
				int map_max = 1;
				map[best_i].submatch_best = emalloc(map_max * sizeof(int));
				if (!map[best_i].submatch_best)
					return 1;
			}
			map[best_i].submatch_best[0] = -1;
			map[best_i].score_with_children = 0;
			map[best_i].compl_z_score = 0;

			for (j = 0; j < best_ctr; j++) {

				if (i == j)
					continue;

				map_complementarity(map + best_i, map + map_best[j], &z);

				map_consistence(NX, NY, &combined_map, map + best_i,
						map + map_best[j], &total_assigned_score, &gap_score);
				consistent = ((map + best_i)->assigned_score
						< total_assigned_score
						&& (map + map_best[j])->assigned_score
								< total_assigned_score);
				if (consistent) {
					score = total_assigned_score;
					if (score > best_score) {
						best_score = score;
						map[best_i].submatch_best[0] = map_best[j];
						map[best_i].score_with_children = total_assigned_score;
						map[best_i].compl_z_score = z;
					}
				}
			}

		}
	}

	free_map(&combined_map);
	return 0;
}
