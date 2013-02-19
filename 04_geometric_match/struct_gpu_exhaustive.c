#include "struct.h"


// GPU version of the function

int find_best_triples_exhaustive_parallel_gpu(Representation* X_rep, Representation* Y_rep, int no_top_rmsd,
        double * best_rmsd, int ** best_triple_x, int ** best_triple_y,
        double **best_quat) {
    // initialization of global array of values
    double ** best_quat_array  = dmatrix(no_top_rmsd * NUM_THREADS, 4);
    int ** best_triple_x_array = intmatrix(no_top_rmsd * NUM_THREADS, 3);
    int ** best_triple_y_array = intmatrix(no_top_rmsd * NUM_THREADS, 3);
    double * best_rmsd_array = (double *) malloc(no_top_rmsd * NUM_THREADS * sizeof (double));

    int cnt;

    for (cnt = 0; cnt < NUM_THREADS * no_top_rmsd; ++cnt) {
        best_rmsd_array[cnt] = BAD_RMSD + 1;
        best_triple_x_array[cnt][0] = -1;
    }

    int i, j, k;
    
    PriorityQueue heap = Initialize(TOP_RMSD);

    int * x_type = X_rep->full_type; // no change
    int NX = X_rep->N_full; // no change
    int NY = Y_rep->N_full;
    int * y_type = Y_rep->full_type;
    //double **x = X_rep->full; // no change
    //double **y = Y_rep->full;
    //int x_triple[3], y_triple[3];
    //double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
    //double rmsd; // 
    //double q_init[4] = {0.0}; // no change
    double ** cmx = X_rep->cm; // no change
    double ** cmy = Y_rep->cm; // no change
    double threshold_dist = THRESHOLD;


    Triples_array x_triple_array = {0}, y_triple_array = {0};
    if (!init_triples_array(&x_triple_array, NX * NX * NX)) exit(1);
    if (!init_triples_array(&y_triple_array, NY * NY * NY)) exit(1);
    
    for (i = 0; i < NX; ++i) {
        for (j = 0; j < NX; ++j) {
            if (i == j) continue;
            if (two_point_distance(cmx[j], cmx[i]) > threshold_dist) continue;
            for (k = 0; k < NX; ++k) {
                if ((k == i) || (k == j)) continue;
                if (two_point_distance(cmx[k], cmx[i]) > threshold_dist) continue;
                if (two_point_distance(cmx[k], cmx[j]) > threshold_dist) continue;

                if (!insert_triple_to_array(x_type, i, j, k, &x_triple_array)) {
                    printf("wrong triple type\n");
                    exit(1);
                }
            }
        }
    }

    for (i = 0; i < NY - 2; ++i) {
        for (j = i + 1; j < NY - 1; ++j) {
            if (two_point_distance(cmy[j], cmy[i]) > threshold_dist) continue;
            for (k = j + 1; k < NY; ++k) {
                if (two_point_distance(cmy[k], cmy[i]) > threshold_dist) continue;
                if (two_point_distance(cmy[k], cmy[j]) > threshold_dist) continue;

                if (!insert_triple_to_array(y_type, i, j, k, &y_triple_array)) {
                    printf("wrong triple type\n");
                    exit(1);
                }
            }
        }
    }


    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.hhh_array, y_triple_array.hhh_array,
			      x_triple_array.hhh_cnt, y_triple_array.hhh_cnt, &heap);
    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.hhs_array, y_triple_array.hhs_array,
			      x_triple_array.hhs_cnt, y_triple_array.hhs_cnt, &heap);
    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.hsh_array, y_triple_array.hsh_array,
			      x_triple_array.hsh_cnt, y_triple_array.hsh_cnt, &heap);
    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.hss_array, y_triple_array.hss_array,
			      x_triple_array.hss_cnt, y_triple_array.hss_cnt, &heap);
    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.shh_array, y_triple_array.shh_array,
			      x_triple_array.shh_cnt, y_triple_array.shh_cnt, &heap);
    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.ssh_array, y_triple_array.ssh_array,
			      x_triple_array.ssh_cnt, y_triple_array.ssh_cnt, &heap);
    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.shs_array, y_triple_array.shs_array,
			      x_triple_array.shs_cnt, y_triple_array.shs_cnt, &heap);
    insert_triple_to_heap_gpu(X_rep, Y_rep, x_triple_array.sss_array, y_triple_array.sss_array,
			      x_triple_array.sss_cnt, y_triple_array.sss_cnt, &heap);
    
    
    free_triples_array(&x_triple_array);
    free_triples_array(&y_triple_array);

    int counter = 0;
    Triple rec;

    while (!IsEmpty(heap)) {
        rec = DeleteMax(heap);
        //printf("%lf\n", rec.rmsd);
        memcpy(*(best_triple_y_array + counter), &rec.triple_x, 3 * sizeof (int));
        memcpy(*(best_triple_x_array + counter), &rec.triple_y, 3 * sizeof (int));
        memcpy(best_rmsd_array  + counter, &rec.rmsd, sizeof (double));
        counter++;
    }

    
    // parallel sort of elements of arrays
    sortTriplets(best_triple_x_array, best_triple_y_array, best_rmsd_array, best_quat_array, no_top_rmsd);

    memcpy(*best_triple_y, *best_triple_y_array, no_top_rmsd * 3 * sizeof (int));
    memcpy(*best_triple_x, *best_triple_x_array, no_top_rmsd * 3 * sizeof (int));
    memcpy(best_rmsd, best_rmsd_array, no_top_rmsd * sizeof (double));


    free_dmatrix(best_quat_array);
    free_imatrix(best_triple_x_array);
    free_imatrix(best_triple_y_array);
    free(best_rmsd_array);
    
    return 0;

}

