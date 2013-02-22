/* 
 * File:   match_triplets.h
 * Author: miles
 *
 * Created on July 20, 2012, 2:19 PM
 */

#ifndef MATCH_TRIPLETS_H
#define	MATCH_TRIPLETS_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "struct.h"
    
# define BAD_RMSD   10.0
# define JACKFRUIT   8
# define NUM_THREADS 8
# define CUTOFF_DNA  3.0

# define MAX_TRIPS  5000

typedef struct {
    int  member[3];
    char fingerprint; /* types of the three vectors and their hand */
} TripleID;

# define TYPE1 1<<3
# define TYPE2 1<<2
# define TYPE3 1<<1
# define HAND  1


# define THRESHOLD 25

    
typedef struct{
    int ** hhh_array; 
    int ** hhs_array;  
    int ** hsh_array; 
    int ** hss_array; 
    int ** shh_array; 
    int ** ssh_array; 
    int ** shs_array; 
    int ** sss_array;
    int hhh_cnt;
    int hhs_cnt;  
    int hsh_cnt;  
    int hss_cnt;  
    int shh_cnt;  
    int ssh_cnt;  
    int shs_cnt;  
    int sss_cnt;     
} Triples_array;

int hand (Representation * X_rep,  int *set_of_directions_x);
int same_hand_triple (Representation * X_rep,  int *set_of_directions_x,
		      Representation * Y_rep, int *set_of_directions_y, int set_size);
int distance_of_nearest_approach (Representation * X_rep, int *set_of_directions_x,
				  Representation * Y_rep, int *set_of_directions_y,
				  int set_size,  double * rmsd_ptr);
int cull_by_dna (Representation *X_rep, int *set_of_directions_x,
		 Representation *Y_rep, int *set_of_directions_y,
		 int set_size, Map *map, double cutoff_rmsd);
    
int insert_triple_to_array(int * triple_types, int triple1, int triple2, int triple3, 
			   Triples_array *triple_array);
    
int insert_triple_to_heap (Representation* X_rep, Representation* Y_rep, int ** x_triple_array,
			  int ** y_triple_array, int x_triple_cnt, int y_triple_cnt, PriorityQueue * heap);
extern int insert_triple_to_heap_gpu(Representation* X_rep, Representation* Y_rep, int ** x_triple_array,
				     int ** y_triple_array, int x_triple_cnt, int y_triple_cnt,
				     PriorityQueue * heap, int no_top_rmsd);

int init_triples_array (Triples_array * triples_array, int size);
int free_triples_array (Triples_array * triples_array);    
int sortTriplets(int ** best_triple_x_array, int ** best_triple_y_array,
		 double * best_rmsd_array, double ** best_quat_array, int top_rmsd);
#ifdef	__cplusplus
}
#endif

#endif	/* MATCH_TRIPLETS_H */

