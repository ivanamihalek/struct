/*
This source code is part of deconSTRUCT,
protein structure database search and backbone alignment application.
Written by Ivana Mihalek, with contributions from Mile Sikic.
Copyright  2012 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.

Contact: ivana.mihalek@gmail.com.
*/

# include "struct.h"

/**************************************************************/
/**************************************************************/
int neighborhood_initialize (Representation *** hood_ptr, int number_of_elements) {

    int i;
    int max_possible_hood_size = number_of_elements -1; // any other members of the representation, but itself
    Representation ** hood = emalloc(number_of_elements*sizeof(Representation *));
    for (i=0; i<number_of_elements; i++) {
	hood[i] = emalloc(sizeof(Representation));
	hood[i]->N_full = 0;
	if ( ! (hood[i]->full  = dmatrix(max_possible_hood_size, 3) )) return 1;
	if ( ! (hood[i]->full_type  = emalloc(max_possible_hood_size*sizeof(int)) )) return 1;
    }

    (*hood_ptr) = hood;
    return 0;
 }
int neighborhood_shutdown (Representation ** hood, int number_of_elements) {

    int i;
    for (i=0; i<number_of_elements; i++) {
	free_dmatrix(hood[i]->full);
	free (hood[i]->full_type );
	free(hood[i]);
    }
    free (hood);
        
    return 0;
 }
/**************************************************************/
/**************************************************************/
int rep_initialize (Representation * rep, Descr * descr  ){

    int N = descr->no_of_elements;
    int i, j;
    int construct_compact (Representation *rep);
    
    rep->N_full = N;
    
    rep->full_no_of_strands = descr->no_of_strands;
    if ( ! (rep->full  = dmatrix(N, 3) )) return 1;
    for (i=0; i<N; i++ ) {	
	for (j=0; j<3; j++) rep->full[i][j] = descr->element[i].p[j];
    }
    if ( ! (rep->cm    = dmatrix(N, 3) )) return 1;
    for (i=0; i<N; i++ ) {
	for (j=0; j<3; j++) rep->cm[i][j]   = descr->element[i].cm[j];
    }
   
    if ( ! (rep->translation  = dmatrix(N,3)) ) return 1;
    if ( ! (rep->transl_norm  = emalloc(N*sizeof(double)))) return 1;

    if ( ! (rep->full_type = emalloc(N*sizeof(int)))) return 1;
    if ( ! (rep->length    = emalloc(N*sizeof(int)))) return 1;
   
    for (i=0; i<N; i++ ) {
	rep->full_type[i] = descr->element[i].type;
	rep->length[i]    = descr->element[i].length;
    }

    if ( ! (rep->compact  = dmatrix (N,3) )) return 1;
    if ( ! (rep->compact_type  = emalloc (N*sizeof (int)) )) return 1;
    
	 
    if ( ! (rep->is_rep_by   = emalloc ( N*sizeof (int))) ) return 1;
    /* 1 for the counter, 1 to represent itself: */
    if ( ! (rep->represents  = intmatrix (N, N )) )return 1;
    if ( construct_compact (rep) ) return 1;

    return 0;
 }

/***********************************************************/
/***********************************************************/
 
 int rep_shutdown  (Representation * rep) {

     free_dmatrix (rep->full);
     free_dmatrix (rep->cm);
    
     free_dmatrix (rep->translation);
     free (rep->transl_norm);

     free (rep->full_type);
     free (rep->length);
     
     free_dmatrix (rep->compact );
     free (rep->compact_type);
     
     free (rep->is_rep_by );

     free_imatrix (rep->represents);
   

     memset (rep, 0, sizeof(Representation) );
     
     return 0;
 }

/******************************************/
int find_cluster_avg (int *cluster, int * strand_id,
		      Representation * rep, int ** nbrs,  double *avg ) {

    int i, j, s;
    int first_cluster_member, cluster_member;
    int cluster_size = cluster[0];
    double norm;

    if (! cluster_size) return 1;
    
    first_cluster_member =  cluster[1]; /* chooses direction */

    for (j=0; j<3; j++) avg[j] = 0.0;
    for ( s=1; s <= cluster_size; s++ ) {
	cluster_member = cluster[s];
	i = strand_id [cluster_member];
	for (j=0; j<3; j++) {
	    avg[j] +=
		nbrs[first_cluster_member][cluster_member]
		*rep->full[i][j];
	}
    }
    norm = 0.0;
    for (j=0; j<3; j++)  {
	avg[j] /=  cluster_size;
	norm   += avg[j]*avg[j];
    }
    norm = sqrt(norm);
    for (j=0; j<3; j++) avg[j] /= norm;


    return 0;

}

/**********************************************************/
/**********************************************************/
/**********************************************************/
/**********************************************************/

int construct_compact (Representation *rep) {

    
     int i, j, new_ctr = 0;

     /******************************************/
     if  (options.merge_cosine > 1.0 ) { /* merging directions not called for */

	 for (i=0; i< rep->N_full; i++) {
	     memcpy (rep->compact[i], rep->full[i], 3*sizeof(double));
	     rep->compact_type[i]   = rep->full_type[i];
	     rep->is_rep_by[i] = i;
	     rep->represents[i][0] = 1; /* this is the count */
	     rep->represents[i][1] = i;
	 }
	 rep->N_compact = rep->N_full;

     } else { 
	 /* cluster directions with the directions such that cos > options.merge_cosine */
  	 int count;
	 int  **cluster;
	 int no_of_clusters, c, s;
	 int cluster_member;
	 int find_cluster_avg (int * cluster, int * strand_id,
		      Representation * rep,  int ** nbrs,  double *avg );
	 int ** element, no_of_elements;   /* strand id's for each sheet */
	 int ** nbrs;  /* neighboring matrix for cluster counting */
	 int  sse_ctr, sse_ctr_2, sse_type, type_ctr[3] = {0};
	 int *cluster_representative, *visited, cluster_id;
	 int *parent_cluster;
	 int total_nbrs, cluster_ctr;
	 double dot, **avg;

	 no_of_elements = rep->N_full;

	 if ( ! (element = intmatrix (3, rep->N_full)))
	     return 1;
	 if ( ! (nbrs    = intmatrix (no_of_elements, no_of_elements) ) )
	     return 1;
	 if ( ! (cluster = intmatrix (no_of_elements+1, no_of_elements+1)))
	     return 1;
	 
	 if ( !( visited  = emalloc ( no_of_elements*sizeof(int))) ) return 1;
	 if ( !( parent_cluster= emalloc (rep->N_full*sizeof(int)))) return 1;
	 if ( !( cluster_representative=emalloc(no_of_elements*sizeof(int))))
	     return 1;
	 
	 if ( ! ( avg = dmatrix ( no_of_elements, 3) ) ) return 1;
	 
	 /*******************************************/
	 /*   extract helices and strands           */
	 /*******************************************/
	 for (sse_ctr=0; sse_ctr < rep->N_full; sse_ctr++ ) {
	     sse_type = rep->full_type[sse_ctr]/2; /* I have the types defined as 1, 2 or 4*/
	     element[sse_type][ type_ctr[sse_type] ] = sse_ctr;
	     type_ctr[sse_type]++;
	     
	     parent_cluster[sse_ctr] = -1;
	 }


	 /***********************************************************/
	 /*   for each type, find nbr matrix and related clusters   */
	 /***********************************************************/

	 cluster_ctr = 0;
	 for ( sse_type = 0; sse_type < 3; sse_type++ ) {

	     /* neighbor matrix */
	     total_nbrs = 0;
	     for (sse_ctr=0; sse_ctr < type_ctr[sse_type]; sse_ctr++) {
		 
		 nbrs [sse_ctr][sse_ctr] = 1;
		 i = element[sse_type][sse_ctr];
		 
		 for (sse_ctr_2=sse_ctr+1; sse_ctr_2< type_ctr[sse_type]; sse_ctr_2++){

		     j = element[sse_type][sse_ctr_2];
		     unnorm_dot (rep->full[i], rep->full[j], &dot);
		     if  (  dot  >= options.merge_cosine ) {
			 nbrs [sse_ctr_2][sse_ctr]= nbrs [sse_ctr][sse_ctr_2] = 1;
			 total_nbrs ++;
		     }
		 }
	     }
	     if (!total_nbrs) continue;
	     
	     /* neighbor matrix */
	     if ( type_ctr[sse_type] == 2 && nbrs [0][1] < 0 ) {
		 continue;
	     } else if (type_ctr[sse_type] == 2 && nbrs [0][1] > 0){
		 no_of_clusters = 1;
		 memset (cluster[0], 0, (no_of_elements+1)*(no_of_elements+1)*sizeof(int) );
		 cluster[1][0] = 2; /* size */
		 cluster[1][1] = 0;
		 cluster[1][2] = 1;
	     } else {
		 /* TODO: does cluster_couter  have a bug for the cluster size == 2? */
		 /* or maybe if there nbrs matrix is all 0 ? */
		 /* maybe it doesn't -   just not sure */ 
		 cluster_counter (type_ctr[sse_type], nbrs, 
				  &no_of_clusters, cluster);
	     }
	     
	     /* for each cluster larger than 1, find a representative */
	     
	     for ( c=1; c <= no_of_elements ; c++) { /* c=0 stores isolated strands */
		 if ( cluster[c][0] < 2 ) continue;
		 
		 for ( s=1; s <= cluster[c][0]; s++ ) {
		     cluster_member = cluster[c][s];
		      /* all sse's in this cluster use the same representative */
		     parent_cluster[cluster_member] = cluster_ctr;
		 }
		 
		 /* find average direction to use as a representative : */
		 find_cluster_avg (cluster[c], element[sse_type],
				   rep, nbrs, avg[cluster_ctr]);
# if 0
		 printf ("type %1d   cluster %2d,  size: %2d\n ",
			sse_type, cluster_ctr, cluster[c][0]);
		 printf ("\tcluster members: " );
		 for ( s=1; s <= cluster[c][0]; s++ ) {
		     cluster_member = cluster[c][s];
		     i = element[sse_type][cluster_member];
		     printf (" %2d:%2d ", cluster_member, i );
		 }
		 printf ("\n");
		 vec_out ( avg[cluster_ctr], 3 , "avg");
		 printf ("\n\n");
# endif
		 cluster_ctr++;
		 
	     } /* end loop over all clusters  */
	     
	 }/* end loop over types (helix or strand) */

	 /*******************************************************/
	 /* assign representatives for each element             */
	 /*******************************************************/
	 new_ctr = 0;
	 cluster_id = -1;
	 memset ( visited, 0, no_of_elements*sizeof(int) );	 
	 
	 for (i=0; i< rep->N_full; i++) {
	     rep->represents[i][0] = 0;
	 }
	 
	 for (i=0; i< rep->N_full; i++) {

	     cluster_id = parent_cluster[i];
	     
	     if (   cluster_id < 0) {
		 
		 memcpy (rep->compact[new_ctr], rep->full[i], 3*sizeof(double));
		 rep->compact_type  [new_ctr] = rep->full_type[i];
		 rep->is_rep_by[i] = new_ctr;
		 new_ctr++;
		 
	     } else if ( ! visited[cluster_id] ) {
		 
		 /* this  cluster we see  for the 1st time */
		 visited[cluster_id] = 1;
		 
		 for (j=0; j<3; j++) rep->compact[new_ctr][j] = avg[cluster_id][j];
		 rep->compact_type  [new_ctr] = rep->full_type[i];

		 rep->is_rep_by[i] = new_ctr;
		 cluster_representative[cluster_id] = rep->is_rep_by[i];
		 new_ctr ++;
		 
	     } else { /* this cluster already seen: keep track of the representative */
		 rep->is_rep_by[i] = cluster_representative[cluster_id];
	     }
	     
	     /* sort the representatives */
	     j = rep->is_rep_by[i];
	     rep->represents[j][0] += 1;
	     count =  rep->represents[j][0];
	     rep->represents[ j ] [count] = i;
	     
	 }

	 rep->N_compact = new_ctr;
	 free_imatrix (element);
	 free_imatrix (nbrs);
	 free_imatrix (cluster);
	 
	 free (visited);
	 free (parent_cluster);
	 free (cluster_representative);
	 
	 free_dmatrix(avg);

     }      
# if 0
     int i_rep, ctr;
     printf ("full set size: %d\n",  rep->N_full);
     for (i=0; i< rep->N_full; i++) {
	 printf ("%2d is represented by ", i);
	 j = rep->is_rep_by[i];
	 printf ("%2d ", j);
	 printf ("\n");
     }
     
     for (i_rep=0; i_rep<rep->N_compact; i_rep++ ) {
	 printf ("%2d  represents ", i_rep);
	 for (ctr=1; ctr <= rep->represents[i_rep][0]; ctr++) {
	     i = rep->represents[i_rep][ctr];
	     printf (" * %2d  *", i);
	 }
	 printf ("\n");
     }
     printf (" effective size: %d \n", rep->N_compact);
     exit (1);
# endif
	 
     return 0;
 }


// serialize Representation and save to a binary file

int rep_save (Representation *rep, char * filename) {
    FILE * fp = fopen(filename, "wb");
    int i,j;
    if (fp == NULL) {
        perror("Failed to open file");
        return EXIT_FAILURE;
    }
    
    int N = rep->N_full;
    fwrite(&rep->N_full, sizeof(int), 1, fp);
    fwrite(&rep->full_no_of_strands, sizeof(int), 1, fp);
    for (i = 0 ; i < N; ++i){
        for (j = 0 ; j < 3; ++j){
        fwrite(&rep->full[i][j], sizeof(double), 1, fp);
        }
    }

    for (i = 0 ; i < N; ++i){
        fwrite(*(rep->cm) + 3*i, sizeof(double), 3, fp);
    }
    
    fwrite(rep->full_type, sizeof(int), N, fp);
    fwrite(rep->length, sizeof(int), N, fp);
    
    fclose(fp);
    return 0;
    
}

// read Representation from a serialized binary file

int rep_read(Representation *rep, char * filename) {
    
    int i, j;
    int N;
    int retv;
    
    FILE * fp = fopen(filename, "rb");
    if (fp == NULL) {
        perror("Failed to open file");
        return EXIT_FAILURE;
    }
    
    
    retv = fread(&N, sizeof(int), 1, fp);
    rep->N_full = N;
    
    retv = fread(&rep->full_no_of_strands, sizeof(int), 1, fp);
    
    if ( ! (rep->full  = dmatrix(N, 3) )) return 1;
    for (i=0; i<N; i++ ) {	
	for (j=0; j<3; j++) {
            retv = fread(&rep->full[i][j], sizeof(double), 1, fp); 
        }
    }
    
    
    if ( ! (rep->cm    = dmatrix(N, 3) )) return 1;
    for (i=0; i<N; i++ ) {
	for (j=0; j<3; j++) {
            retv = fread(&rep->cm[i][j], sizeof(double), 1, fp);
        }
    }
   
    if ( ! (rep->translation  = dmatrix(N,3)) ) return 1;
    if ( ! (rep->transl_norm  = emalloc(N*sizeof(double)))) return 1;

    if ( ! (rep->full_type = emalloc(N*sizeof(int)))) return 1;
    if ( ! (rep->length    = emalloc(N*sizeof(int)))) return 1;
   
    retv = fread(rep->full_type, sizeof(int), N, fp);
    retv = fread(rep->length, sizeof(int), N, fp);

    if ( ! (rep->compact  = dmatrix (N,3) )) return 1;
    if ( ! (rep->compact_type  = emalloc (N*sizeof (int)) )) return 1;
    
	 
    if ( ! (rep->is_rep_by   = emalloc ( N*sizeof (int))) ) return 1;
    /* 1 for the counter, 1 to represent itself: */
    if ( ! (rep->represents  = intmatrix (N, N )) )return 1;
    if ( construct_compact (rep) ) return 1;

    fclose(fp);
    return 0;
    
    
}

