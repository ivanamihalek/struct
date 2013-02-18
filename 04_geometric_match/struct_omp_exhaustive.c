// old OpenMP version

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
# ifdef OMP
int find_best_triples_exhaustive_parallel(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
        double **best_quat) {
    // initialization of global array of values
    no_top_rmsd = TOP_RMSD;
    // printf("%d\n", no_top_rmsd);

    // printf("proba %d %lf\n", X_rep->N_full, X_rep->cm[0][0]);

    double ** best_quat_array = dmatrix(no_top_rmsd * NUM_THREADS, 4);
    int ** best_triple_x_array = intmatrix(no_top_rmsd * NUM_THREADS, 3);
    int ** best_triple_y_array = intmatrix(no_top_rmsd * NUM_THREADS, 3);
    double * best_rmsd_array = (double *) malloc(no_top_rmsd * NUM_THREADS * sizeof (double));

    int cnt;

    for (cnt = 0; cnt < NUM_THREADS * no_top_rmsd; ++cnt) {
        best_rmsd_array[cnt] = BAD_RMSD + 1;
        best_triple_x_array[cnt][0] = -1;
    }

    omp_set_num_threads(NUM_THREADS);

#pragma omp parallel 

    {

        int top_ctr, i, j, k, l, n, m;
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
        double threshold_dist = THRESHOLD;

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
        for (i = 0; i < NX; ++i) {
            for (j = 0; j < NY - 2; ++j) {
                if (x_type[i] != y_type[j]) continue;
                
                for (k = 0; k < NX; ++k) {
                    if (k == i) continue;
                    if (two_point_distance(cmx[i], cmx[k]) > THRESHOLD) continue;
                    
                    for (l = j + 1; l < NY - 1; ++l) {
                        if (x_type[k] != y_type[l]) continue; 
                        if (two_point_distance(cmy[j], cmy[l]) > THRESHOLD) continue;
                            
                        for (m = 0; m < NX; ++m) {
                            if (m == k || m == i) continue;
                                
                            if (two_point_distance(cmx[i], cmx[m]) > THRESHOLD) continue;
                            if (two_point_distance(cmx[k], cmx[m]) > THRESHOLD) continue;
                            
                            for (n = l + 1; n < NY; ++n) {
                                if (x_type[m] != y_type[n]) continue;
                                if (two_point_distance(cmy[j], cmy[n]) > THRESHOLD) continue;
                                if (two_point_distance(cmy[l], cmy[n]) > THRESHOLD) continue;
                                    
                                x_triple[0] = i;
                                y_triple[0] = j;
                                x_triple[1] = k;
                                y_triple[1] = l;
                                x_triple[2] = m;
                                y_triple[2] = n;

                                
                                
                                if (!same_hand_triple(X_rep, x_triple, Y_rep, y_triple, 3)) continue;
                                if (distance_of_nearest_approach(X_rep, x_triple,
                                        Y_rep, y_triple, 3, &rmsd))     continue;
                                if (rmsd > cutoff_rmsd)     continue;

                                //if (opt_quat(x, NX, x_triple, y, NY, y_triple, 3, q_init, &rmsd)) continue;

                                
                                
                                for (top_ctr = 0; top_ctr < no_top_rmsd; top_ctr++) {
                                    // insertion of a new values in arrays keeping arrays sorted

                                    if (rmsd <= best_rmsd_local[top_ctr]) {
                                        chunk = no_top_rmsd - top_ctr - 1;

                                        if (chunk) {
                                            memmove(best_rmsd_local + top_ctr + 1,
                                                    best_rmsd_local + top_ctr, chunk * sizeof(double));

                                            memmove(best_quat_local[top_ctr + 1],
                                                    best_quat_local[top_ctr], chunk * 4 * sizeof(double));

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

            //           printf("%d\n", i);

            // each thread copies values to global arrays in accordance with its thread id
            //           printf("myid: %d\n", myid);

            memcpy(*(best_quat_array + myid * no_top_rmsd), *(best_quat_local), no_top_rmsd * 4 * sizeof (double));

            memcpy(*(best_triple_y_array + myid * no_top_rmsd), *(best_triple_y_local), no_top_rmsd * 3 * sizeof (int));
            memcpy(*(best_triple_x_array + myid * no_top_rmsd), *(best_triple_x_local), no_top_rmsd * 3 * sizeof (int));
            memcpy(best_rmsd_array + myid*no_top_rmsd, best_rmsd_local, no_top_rmsd * sizeof (double));

        }

        free_dmatrix(best_quat_local);
        free_imatrix(best_triple_x_local);
        free_imatrix(best_triple_y_local);
        free(best_rmsd_local);

        // parallel sort of elements of arrays 
        sortTriplets(best_triple_x_array, best_triple_y_array, best_rmsd_array, best_quat_array, no_top_rmsd);

    }


    // 
    /*
        memcpy(*best_quat, *best_quat_array, no_top_rmsd * 4 * sizeof(double));
     */
    memcpy(*best_triple_y, *best_triple_y_array, no_top_rmsd * 3 * sizeof (int));
    memcpy(*best_triple_x, *best_triple_x_array, no_top_rmsd * 3 * sizeof (int));
    memcpy(best_rmsd, best_rmsd_array, no_top_rmsd * sizeof (double));


    free_dmatrix(best_quat_array);
    free_imatrix(best_triple_x_array);
    free_imatrix(best_triple_y_array);
    free(best_rmsd_array);

    return 0;

}
# endif
