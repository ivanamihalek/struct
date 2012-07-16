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

/**************************************************************/
/*  postprocess - find the actual tf                          */
/* 	   & mapping on the bb level                          */
/**************************************************************/

/**************************************************************/
int following_loop (int *element_begin, int *element_end,
		    int no_of_elements, int no_of_res, 
		    int element_ctr, int * first_res, int * last_res);

int preceding_loop (int *element_begin, int *element_end,
		    int element_ctr, int * first_res, int * last_res);


int single_map_align_backbone (Descr *descr1, Protein * protein1, Representation *rep1, 
		    Descr *descr2, Protein * protein2, Representation *rep2, 
		    Map *map);


int align_backbone (Descr *descr1, Protein * protein1, Representation *rep1, 
		    Descr *descr2, Protein * protein2, Representation *rep2, 
		    List_of_maps *list){
    
     if ( list->best_array_used == 0) return 1;
     
     int best_ctr, map_ctr, retval;
     int * new_best;
     double *bb_score_array;
     Map *current_map;

     if ( !(bb_score_array = emalloc(list->best_array_used*sizeof(double))))  return 1;    
     if ( !(new_best = emalloc(list->best_array_used*sizeof(int))))  return 1;    
 
     
     for (best_ctr=0; best_ctr<list->best_array_used; best_ctr++) {

	 map_ctr = list->map_best[best_ctr];
	 if ( map_ctr >= list->best_array_allocated) {
	     printf ("%s:%d error assigning array size.\n", __FILE__, __LINE__ );
	     exit (1);
	 }
	 current_map = list->map+map_ctr;
	 retval = single_map_align_backbone (descr1, protein1, rep1, descr2, protein2, rep2, current_map);
	 if (retval) {
	     printf (" error doing bb alignment   db:%s  query:%s \n",
		     descr1->name, descr2->name);
	     exit (retval);
	 }
	 
	 new_best[best_ctr] = map_ctr;
	 /* the array_qsort sorts in the increasing order, so use negative aln score: */
	 bb_score_array[best_ctr] = -current_map->aln_score;
     }

     array_qsort ( new_best, bb_score_array, list->best_array_used);

     memcpy (list->map_best, new_best, list->best_array_used*sizeof(int));
     
     free (bb_score_array);
     free (new_best);
     
     return 0;
}

/**************************************************************/
int single_map_align_backbone (Descr *descr1, Protein * protein1, Representation *rep1, 
			       Descr *descr2, Protein * protein2, Representation *rep2, 
			       Map * map) {
    
    /* for now,  we will just postprocess the  best map */

    int no_res_1 = protein1->length, no_res_2= protein2->length;
    int resctr1, resctr2;
    int map_size;
    int *residue_map_i2j, *residue_map_j2i;
    
    int *type_1, *type_2;
    int longest_element_length = (protein1->length > protein2->length) ?
	    protein1->length : protein2->length;
   
    double d0 = options.distance_tol_in_bb_almt;
    double aln_score, rmsd;
    double ** similarity;
    double ** sim_in_element;
    double **x, **y;
    double **R, T[3], q[4];
    double total_score = 0;
     /* for the MC: */
    int max_no_steps = 20, no_steps = 0;
    int done = 0, toggle = 0;
    double *current_q,  *old_q, *current_T, *old_T;
    double *best_q, *best_T;
    double old_score = total_score, current_score = 0.0, d_mc = 0.5;
    double max_score;
    double t_mc, d_init;

      
    int  alignment_size  (int * residue_map_i2j, int no_res_1 );
    int  closeness_score_for_sse_almt (Descr *descr1, Representation *rep1, Representation *rep2, Map * map,
			  Protein *protein1, Protein *protein2,
			  double **R, double *T, double d0, double ** similarity, double * score_ptr);
    int following_loop (int *element_begin, int *element_end,
			int no_of_elements, int no_of_res, 
			int element_ctr, int * first_res, int * last_res);  
    int map2rotation (Protein *protein1, Protein *protein2, int *residue_map_i2j,
		       double **x, double **y, double *q, double *T, double *rmsd);
    int out_of_order_alignment (Descr *descr1,  Descr *descr2, Map *map, int *element_1_begin, int *element_1_end,
			     int *element_2_begin, int *element_2_end,
			     int longest_element_length, double ** similarity, double ** sim_in_element,
				 int *residue_map_i2j, int *residue_map_j2i, double * score_ptr);
    int preceding_loop (int *element_begin, int *element_end,
			int element_ctr, int * first_res, int * last_res);
    
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    
    construct_translation_vecs (rep1, rep2, map);
    
    /* make sure that we have all the info we might need */
  
    /* define matrix, the size of nr of residues in set of SSEs 
       x nr of  residues in the other set of SSEs, and fill it with -1 */
    similarity = dmatrix (no_res_1, no_res_2);
    if ( !similarity ) return 1;
    
    sim_in_element = dmatrix (no_res_1, no_res_2);
    if ( !similarity ) return 1;
    
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	for (resctr2=0; resctr2<no_res_2; resctr2++) {
	    similarity[resctr1][resctr2] = -1;
	}
    }
    
    /* alloc */

    type_1 = protein1->sse_sequence;
    type_2 = protein2->sse_sequence;
    
    
    if ( ! (residue_map_i2j     = emalloc (no_res_1*sizeof(int))) ) return 1;
    if ( ! (residue_map_j2i     = emalloc (no_res_2*sizeof(int))) )  return 2;

    if ( ! (x = dmatrix (3, no_res_1+no_res_2)))  exit(1);
    if ( ! (y = dmatrix (3, no_res_1+no_res_2)))  exit(1);

    
    /*********************************************************************/
    /*********************************************************************/
    /* aliases */
    int *element_1_begin, *element_1_end; /* "element" here means SSE */
    int *element_2_begin, *element_2_end;
 
    element_1_begin = protein1->element_begin;
    element_1_end   = protein1->element_end;
     
    element_2_begin = protein2->element_begin;
    element_2_end   = protein2->element_end;
 


    /************************************************************/
    /************************************************************/
    /* ALIGNMENT, round 1                                       */
    /************************************************************/
    /* for all mapped blocks calculate similarity as exp (-d/d0) */
    total_score = 0.0;
    quat_to_R (map->q, R);

    closeness_score_for_sse_almt (descr1, rep1, rep2, map, protein1, protein2,
		     R, NULL, d0, similarity, &total_score);

    /* run Smith-Waterman and use the mapped CA to find the transformation        */
    /* I have another copy of SW here (the first one is in struct_map.c)          */
    /* so I wouldn't fumble with parameters - the two should be joined eventually */
    if ( options.search_algorithm == SEQUENTIAL ) {

	smith_waterman_2 (no_res_1, no_res_2, similarity,
			  residue_map_i2j, residue_map_j2i, &aln_score);
    } else {
	
	out_of_order_alignment (descr1, descr2, map, element_1_begin, element_1_end,
				element_2_begin, element_2_end, longest_element_length,
				similarity, sim_in_element,
				residue_map_i2j,  residue_map_j2i, &aln_score);
    }
 
    
	
    map2rotation (protein1, protein2, residue_map_i2j, x, y, q, T, &rmsd);
    
    quat_to_R (q, R);
    current_score = alignment_score (protein1, protein2, residue_map_i2j, R, T, d0);

    /*********************************************************/
    /* fiddle iteratively with the transformation            */
    if ( ! (current_q = emalloc (4*sizeof(double)) )) return 1;
    if ( ! (old_q     = emalloc (4*sizeof(double)) )) return 1;
    if ( ! (best_q    = emalloc (4*sizeof(double)) )) return 1;
    if ( ! (current_T = emalloc (3*sizeof(double)) )) return 1;
    if ( ! (old_T     = emalloc (3*sizeof(double)) )) return 1;
    if ( ! (best_T    = emalloc (3*sizeof(double)) )) return 1;
    
    srand48 (time (0));
    
    memcpy (current_q, q, 4*sizeof(double));
    memcpy (    old_q, q, 4*sizeof(double));
    memcpy (   best_q, q, 4*sizeof(double));
    
    memcpy (    old_T, T, 3*sizeof(double));
    memcpy (   best_T, T, 3*sizeof(double));
    memcpy (current_T, T, 3*sizeof(double));

    quat_to_R ( current_q, R);

    d_init = d0;


    /* t_mc = exp ( (1.0- (double)anneal_round)/10.0); */
    d_mc = d_init;
    t_mc = 5;
	
    memcpy (current_q, best_q, 4*sizeof(double));
    memcpy (current_T, best_T, 3*sizeof(double));
    memcpy (old_q, best_q, 4*sizeof(double));
    memcpy (old_T, best_T, 3*sizeof(double));
    old_score = 0;
    max_score = 0;
    no_steps  = 0;
    
    toggle = 1;
    done = 0;
    
    while (no_steps < max_no_steps && !done ) {

	    
	closeness_score_for_sse_almt (descr1, NULL, NULL, map, protein1, protein2,
			 R, current_T, d_mc, similarity, &total_score);
	    
	
	if ( options.search_algorithm == SEQUENTIAL ) {
	    smith_waterman_2 (no_res_1, no_res_2, similarity,
			      residue_map_i2j,  residue_map_j2i, &aln_score);
	} else {
	    out_of_order_alignment (descr1, descr2, map, element_1_begin, element_1_end,
				    element_2_begin,  element_2_end, longest_element_length,
				    similarity, sim_in_element,
				    residue_map_i2j,  residue_map_j2i, &aln_score);
	}
	
 	map2rotation (protein1, protein2, residue_map_i2j, x, y, current_q,  current_T, &rmsd);
	
	quat_to_R ( current_q, R);
	current_score = alignment_score (protein1, protein2, residue_map_i2j, R, current_T, d_mc);

	
	if ( current_score >  max_score )  {
	    max_score = current_score;
	    memcpy (best_q, current_q, 4*sizeof(double));
	    memcpy (best_T, current_T, 3*sizeof(double));
	}

	if (old_score) done = ( fabs(old_score-current_score)/old_score < 0.01);
	old_score = current_score;
	no_steps++;
    
    }

    
    memcpy (q, best_q, 4*sizeof(double));
    memcpy (T, best_T, 3*sizeof(double));
	 
    

    free (current_q);
    free (old_q);
    free (best_q);
    free (current_T);
    free (old_T);
    free (best_T);
   


    /************************************************************/
    /************************************************************/
    /* ALIGNMENT, round 2                                       */
    /************************************************************/
    /************************************************************/
    /* find the similarity matrix for this new rotation  -- this*/
    /* time extending to neighboring elements                   */

    closeness_score_for_bb_almt (map, protein1, protein2, R, T, d0,
				 similarity, &total_score);

    /************************************************************/

    
    memset (residue_map_i2j, 0, no_res_1*sizeof(int)); 
    memset (residue_map_j2i, 0, no_res_2*sizeof(int)); 

    if ( options.search_algorithm == SEQUENTIAL ) {
	smith_waterman_2 (no_res_1, no_res_2, similarity,
			  residue_map_i2j,  residue_map_j2i, &aln_score);
    } else {
	out_of_order_alignment (descr1, descr2, map, element_1_begin, element_1_end,
				element_2_begin,  element_2_end, longest_element_length,
				similarity, sim_in_element,
				residue_map_i2j,  residue_map_j2i, &aln_score);
    }

    map2rotation (protein1, protein2, residue_map_i2j, x, y, q, T, &rmsd);
    
    quat_to_R (q, R);
    
    //aln_score = alignment_score (protein1, protein2, residue_map_i2j, R, T, d0);
    map_size  = alignment_size  (residue_map_i2j, protein1->length);		

    memcpy (&(map->q[0]), &q[0], 4*sizeof(double));
    memcpy (&(map->T[0]), &T[0], 3*sizeof(double));
    
    map->x2y_residue_level  = residue_map_i2j;
    map->y2x_residue_level  = residue_map_j2i;
    
    map->x2y_residue_l_size = no_res_1;
    map->y2x_residue_l_size = no_res_2;
    

    /*************************************************************************/
    map->res_almt_length     = map_size;
    map->aln_score           = aln_score;
    map->res_rmsd            = rmsd;

    free_dmatrix(R);
    free_dmatrix (similarity);
    free_dmatrix (sim_in_element);
    free_dmatrix (x);
    free_dmatrix (y);

 
     
    return 0;
}


/*************************************************************************/
/*************************************************************************/
int  alignment_size ( int * residue_map_i2j, int no_res_1 ) {

    int resctr1, resctr2;
    int aln_size =0;
    
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = residue_map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	aln_size++;
    }
    return aln_size;
}

/*************************************************************************/
/*************************************************************************/
double alignment_score (Protein * protein1, Protein * protein2, int * residue_map_i2j,
			double **R, double *T,  double d0) {

    double d;
    double aln_score = 0.0;
    double ca1[3], ca2[3], rotated_ca1[3];
    int resctr1, resctr2;
    int no_res_1 = protein1->length;
    int aln_size =0;
    
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = residue_map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	find_Calpha ( protein1, resctr1, ca1 );
	find_Calpha ( protein2, resctr2, ca2 );
	
	/* rotate and  translate*/
	point_rot_tr (ca1, R, T, rotated_ca1);
		
	d = two_point_distance (rotated_ca1, ca2);
	aln_score +=  exp (-d/d0);

	aln_size++;
    }
    return aln_score;
}

/*************************************************************************/
/*************************************************************************/
/* the sse-s are matched by their spatial position, and a rough alignment
   exists on the bb level */
int  closeness_score_for_bb_almt (Map *map,  Protein *protein1, Protein *protein2,
				  double **R, double *T, double d0, double **similarity, double *score_ptr)
{
    int element_ctr_1, element_ctr_2;
    int resctr1, resctr2;
    double d, score;
    int *element_1_begin   = protein1->element_begin;
    int *element_1_end     = protein1->element_end;
    int *element_2_begin   = protein2->element_begin;
    int *element_2_end     = protein2->element_end;
    int no_res_1 = protein1->length;
    int no_res_2 = protein2->length;
    

    int last_res_prev_loop_1, last_res_prev_loop_2; 
    int first_res_prev_loop_1, first_res_prev_loop_2;
    int last_res_1, last_res_2; 
    int first_res_1, first_res_2;
    int last_res_next_loop_1, last_res_next_loop_2; 
    int first_res_next_loop_1, first_res_next_loop_2;
    int last_element_1, last_element_2;

    int no_of_elements_1 = protein1->no_helices + protein1->no_strands;
    int no_of_elements_2 = protein2->no_helices + protein2->no_strands;
    
    double ca1[3], ca2[3], rotated_ca1[3];
    double aux; // step = NR_POINTS/MAX_EXP_VALUE;
   /************/
    
    for (resctr1=0; resctr1<protein1->length; resctr1++) {
	for (resctr2=0; resctr2<protein2->length; resctr2++) {
	    similarity[resctr1][resctr2] = -1;
	}
    }
  
    /* for all mapped blocks calculate similarity as exp (-d/d0) */
    score = 0.0;
    
    
    last_element_1 = 0;
    last_element_2 = 0;
    for (element_ctr_1=0; element_ctr_1 <  no_of_elements_1; element_ctr_1++) {
	
	element_ctr_2 = map->x2y[element_ctr_1];
	if (element_ctr_2 <  0) continue;

	last_element_1 = element_ctr_1;
	last_element_2 = element_ctr_2;
	

	/****************************************************/
	/* try to extend the match to the previous loop     */

	/* the first and the last residues of the preceding loops */
	if ( element_1_begin[element_ctr_1] == 0 ) {
	     first_res_prev_loop_1 = 0;
	} else {
	     preceding_loop (element_1_begin, element_1_end, element_ctr_1,
			     &first_res_prev_loop_1, &last_res_prev_loop_1);
	}
	
	if ( element_2_begin[element_ctr_2] == 0 ) {
	     first_res_prev_loop_2 = 0;
	} else {
	     preceding_loop (element_2_begin, element_2_end, element_ctr_2,
			     &first_res_prev_loop_2, &last_res_prev_loop_2);
	}

	/****************************************************/
	/* now try to extend the match to the following loop */
	/* the first and the last residues of the following loops */
	/* there is some redundancy here (prev/next) but we can */
        /* get rid of it when it becomes a bottleneck */
	if ( element_1_end[element_ctr_1] >= no_res_1 -1) {
	     last_res_next_loop_1 = no_res_1 -1;
	} else {
	     following_loop (element_1_begin, element_1_end, no_of_elements_1, no_res_1,
			     element_ctr_1, &first_res_next_loop_1, &last_res_next_loop_1);
	}
	
	if ( element_2_end[element_ctr_2] >= no_res_2 -1) {
	     last_res_next_loop_2 = no_res_2 -1;
	} else {
	     following_loop (element_2_begin, element_2_end, no_of_elements_2, no_res_2,
			     element_ctr_2, &first_res_next_loop_2, &last_res_next_loop_2);
	}

	/* this is an exact map point  -- these two we won't try aligning to anything else */
	for (resctr1 = first_res_prev_loop_1; resctr1<= last_res_next_loop_1; resctr1++) {
	    /* find the representative atom (CA) for residue resctr1*/
	    if ( find_Calpha ( protein1, resctr1, ca1 ) ) {
		/* find_Calpha returns 1 on failure */
		for (resctr2= first_res_prev_loop_2; resctr2<= last_res_next_loop_2; resctr2++)
		    similarity[resctr1][resctr2] = 0.0;
		continue;
	    }
	    /* rotate and  translate*/
	    point_rot_tr (ca1, R, T, rotated_ca1);
		
	    for (resctr2= first_res_prev_loop_2; resctr2<= last_res_next_loop_2; resctr2++) {
		/* skip if we are sure the two are not of the same type */
		if ( (protein1->sse_sequence[resctr1]|protein2->sse_sequence[resctr2])
		     == (HELIX|STRAND) )  continue;
	

		/* find the representative atom (CA) for residue resctr2*/
		if ( find_Calpha (protein2, resctr2,  ca2 ) ) {
		    similarity[resctr1][resctr2] = 0.0;
		    continue;
		}
		d = two_point_distance (rotated_ca1, ca2);
		
		if ( d<options.max_almt_dist) { 
		    aux = d/d0;
		    if (  aux <  MAX_EXP_VALUE ) {
			//similarity[resctr1][resctr2] = exp_table[(int)(aux*step)];
			similarity[resctr1][resctr2] = exp (-aux);
			score += similarity[resctr1][resctr2];
		    } else {
			similarity[resctr1][resctr2] = 0.0;
		    }
		}
	    }
	}



    }
    
    /* another ad-hoc fix: did we skip an element of the
       same type in both sequences? */
    int gap = 0;
    int gap_begin_1 = 0, gap_end_1 = 0;
    int gap_begin_2 = 0, gap_end_2 = 0;
    int all_healthy;
    for (element_ctr_1=0; element_ctr_1 < no_of_elements_1; element_ctr_1++) {
     	element_ctr_2 = map->x2y[element_ctr_1];

	if (element_ctr_2>=0) {
	  if (gap && gap <=2) {
	    gap_end_1 = element_ctr_1-1;
	    gap_end_2 = element_ctr_2-1;

	    all_healthy  = (gap_begin_1 >=0);
	    all_healthy &= (gap_end_1 >=0);
	    all_healthy &= (gap_begin_2 >=0);
	    all_healthy &= (gap_end_2 >=0);
	    all_healthy &= ( gap_end_1 - gap_begin_1 <= 2);
	    all_healthy &= ( gap_end_2 - gap_begin_2 <= 2);

	    if ( all_healthy ) {
		first_res_1  = element_1_begin[gap_begin_1];
		last_res_1   = element_1_end  [gap_end_1];
		first_res_2  = element_2_begin[gap_begin_2];
		last_res_2   = element_2_end  [gap_end_2];
		for (resctr1=first_res_1; resctr1<=last_res_1; resctr1++) {
		    if ( find_Calpha ( protein1, resctr1, ca1 ) ) {
			/* find_Calpha returns 1 on failure */
			for (resctr2= first_res_prev_loop_2; resctr2<= last_res_next_loop_2; resctr2++)
			    similarity[resctr1][resctr2] = 0.0;
			continue;
		    }
		    /* rotate and  translate*/
		    point_rot_tr (ca1, R, T, rotated_ca1);
		    /* the two, however, cannot be matched if they are not
		       of the same type (helix, strand or unassigned)*/
		    for (resctr2=first_res_2; resctr2<=last_res_2; resctr2++) {
			
			if ( (protein1->sse_sequence[resctr1]|
			      protein2->sse_sequence[resctr2])
			     == (HELIX|STRAND) )  continue;
			
			if ( find_Calpha (protein2, resctr2,  ca2 ) ) {
			    similarity[resctr1][resctr2] = 0.0;
			    continue;
			}

			d = two_point_distance (rotated_ca1, ca2);
			if ( d< options.max_almt_dist) { 
			    aux = d/d0;
			    if (  aux <  MAX_EXP_VALUE ) {
				//similarity[resctr1][resctr2] = exp_table[(int)(aux*step)];
				similarity[resctr1][resctr2] = exp (-aux);
				score += similarity[resctr1][resctr2];
			    } else {
				similarity[resctr1][resctr2] = 0.0;
			    }
			}

		    }
		}
	    }
	    gap = 0;
	    gap_begin_1 = element_ctr_1+1;
	    gap_begin_2 = element_ctr_2+1;
	    continue;
	  }
	  gap++;
	}
    }
    



    
    *score_ptr = score; 
    return 0;

    
}


/*************************************************************************/
/*************************************************************************/
/* the sse-s are matched by their spatial position, and a rough alignment
   exists on the bb level */
int  closeness_score_for_bb_almt_within_sse_only (Map *map,  Protein *protein1, Protein *protein2,
				  double **R, double *T, double d0, double **similarity, double *score_ptr)
{
    int element_ctr_1, element_ctr_2;
    int resctr1, resctr2;
    double d, score;
    int *element_1_begin   = protein1->element_begin;
    int *element_1_end     = protein1->element_end;
    int *element_2_begin   = protein2->element_begin;
    int *element_2_end     = protein2->element_end;

    int no_of_elements_1 = protein1->no_helices + protein1->no_strands;
    
    double ca1[3], ca2[3], rotated_ca1[3];
    double aux; // step = NR_POINTS/MAX_EXP_VALUE;
    /************/
    
    for (resctr1=0; resctr1<protein1->length; resctr1++) {
	for (resctr2=0; resctr2<protein2->length; resctr2++) {
	    similarity[resctr1][resctr2] = -1;
	}
    }
  
    /* for all mapped blocks calculate similarity as exp (-d/d0) */
    score = 0.0;
    
    for (element_ctr_1=0; element_ctr_1 < no_of_elements_1; element_ctr_1++) {

	/* the corresponding element in the other structure: */
	element_ctr_2 = map->x2y[element_ctr_1];
	if (element_ctr_2 < 0) continue;
	 
	for (resctr1=element_1_begin[element_ctr_1]; 
	     resctr1<= element_1_end[element_ctr_1]; resctr1++) {
	    
	    /* find the representative atom (CA) for residue resctr1*/
	    find_Calpha ( protein1, resctr1, ca1 );
	    
	    point_rot_tr (ca1, R, T, rotated_ca1);
	    
	    for (resctr2=element_2_begin[element_ctr_2]; 
		 resctr2<= element_2_end[element_ctr_2]; resctr2++) {
		  
		/* find the representative atom (CA) for residue resctr2*/
		find_Calpha (protein2, resctr2,  ca2);
		   
		/* finally, find the distance & assign the "similarity" score */
		d = two_point_distance (rotated_ca1, ca2);
		
		if ( d< options.max_almt_dist ) { 
		    aux = d/d0;
		    if (  aux <  MAX_EXP_VALUE ) {
			//similarity[resctr1][resctr2] = exp_table[(int)(aux*step)];
			similarity[resctr1][resctr2] = exp (-aux);
		       } else {
			   similarity[resctr1][resctr2] = 0.0;
		       }
		   }
		   score += similarity[resctr1][resctr2];
	      }
	 }
    }

    *score_ptr = score; 
    return 0;

    
}








/*************************************************************************/
/*************************************************************************/
/* the relative positions of sse-s need to be reconstituted in space      */
int  closeness_score_for_sse_almt (Descr *descr1,
		      Representation *rep1, Representation *rep2, Map *map, 
		      Protein *protein1, Protein *protein2,
		      double **R, double *fixed_T, double d0, double ** similarity, double * score_ptr) {

    int element_ctr_1, element_ctr_2;
    int resctr1, resctr2, i, j;
    double d, score;
    int *element_1_begin   = protein1->element_begin;
    int *element_1_end     = protein1->element_end;
    int *element_2_begin   = protein2->element_begin;
    int *element_2_end     = protein2->element_end;
    
    double T[3];
    double ca1[3], ca2[3], rotated_ca1[3];
    double aux; // step = NR_POINTS/MAX_EXP_VALUE;
    /************/
    int find_Calpha ( Protein *protein, int  resctr, double ca[3] );
    
    if ( !rep1 )  {
	if ( !fixed_T) {
	    fprintf (stderr, "Error in closeness_score_for_sse_almt().\n");
	    exit (1);
	}
	memcpy (T, fixed_T, 3*sizeof(double));
    } else if (!rep2 ) {
	fprintf (stderr, "Error in closeness_score_for_sse_almt() -- both reps or none ...\n");
	exit (1);
    }
  
    /* for all mapped blocks calculate similarity as exp (-d/d0) */
    score = 0.0;
     
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {
	
	element_ctr_2 = map->x2y[element_ctr_1];
	if (element_ctr_2 < 0) continue;

	if ( rep1 ) { /* "exploding" from the origin */
	    for (i=0; i<3; i++) {
		T[i] = 0.0;
		for (j=0; j<3; j++) T[i] += R[i][j]*rep1->translation[element_ctr_1][j];
	    }
	}

	
	for (resctr1=element_1_begin[element_ctr_1]; 
	     resctr1<= element_1_end[element_ctr_1]; resctr1++) {
	    
	    /* find the representative atom (CA) for residue resctr1*/
	    find_Calpha ( protein1, resctr1, ca1 );
	    

	    if ( rep1) {/* translate to the origin */
		for (i=0; i<3; i++) ca1[i] -= rep1->origin[i];
	    }
	      
	    /* rotate and  translate out*/
	    point_rot_tr (ca1, R, T, rotated_ca1);
	    
	    for (resctr2=element_2_begin[element_ctr_2]; 
		 resctr2<= element_2_end[element_ctr_2]; resctr2++) {
		  
		/* find the representative atom (CA) for residue resctr2*/
		find_Calpha (protein2, resctr2,  ca2);
		   
		/* translate to cm & out */
		if ( rep2 ) {
		    for (i=0; i<3; i++) ca2[i] += -rep2->origin[i] + rep2->translation[element_ctr_2][i];
		}
		/* finally, find the distance & assign the "similarity" score */
		d = two_point_distance (rotated_ca1, ca2);		   
		if ( d<options.max_almt_dist) { 
		    aux = d/d0;
		    if (  aux <  MAX_EXP_VALUE ) {
			//similarity[resctr1][resctr2] = exp_table[(int)(aux*step)];
			similarity[resctr1][resctr2] = exp (-aux);
		    } else {
			similarity[resctr1][resctr2] = 0.0;
		    }
		}

		score += similarity[resctr1][resctr2];
	    }
	}

    }
    
    *score_ptr = score;
    
    return 0;
}


/*************************************************************************/
/*************************************************************************/
int preceding_loop (int *element_begin, int *element_end,
		    int element_ctr, int * first_res, int * last_res) {
    
    if ( element_ctr == 0 ) {
	/* there is no previous element */
	*first_res = 0;
    } else {
	*first_res = element_end[element_ctr-1]+1;
    }
    *last_res = element_begin[element_ctr]-1;

    return 0;
}

/*************************************************************************/
/*************************************************************************/
int following_loop (int *element_begin, int *element_end,
		    int no_of_elements, int no_of_res, 
		    int element_ctr, int * first_res, int * last_res) {
    
    *first_res = element_end[element_ctr]+1;
    if (element_ctr >= no_of_elements-1) {
	*last_res = no_of_res - 1;
    } else {
	*last_res = element_begin[element_ctr+1]-1;
    }

    return 0;
}

/*************************************************************************/
/*************************************************************************/
double  two_point_distance (double point1[3], double point2[3] ) {
    int i;
    double aux;
    double d = 0;
    for (i=0; i<3; i++) {
	aux = point1[i] - point2[i];
	d  += aux*aux;
    }
    return  sqrt (d);
}
/*************************************************************************/
/*************************************************************************/
int point_rot_tr (double point_in[3], double **R, double T[3],double point_out[3]) {

    int i, j;
    /* rotate */
    for (i=0; i<3; i++) {
	point_out[i] = 0.0;
	for (j=0; j<3; j++) point_out[i] += R[i][j]*point_in[j];
    }
    /* translate out */
    for (i=0; i<3; i++) point_out[i] += T[i];
    return 0;
}
/*************************************************************************/
/*************************************************************************/
int smith_waterman_2 (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char   **direction;
    double gap_opening   = -0.2;
    double gap_extension = -0.1;
    double endgap        =  0.0;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    double F_max;
    int use_endgap = 1;
    int F_max_i, F_max_j;
    int i,j;

    /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;


    
    /* fill the table */
    F_max = FAR_FAR_AWAY;
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
    for( i= max_i; i>F_max_i; i--) map_i2j[i-1] = FAR_FAR_AWAY;;
    for( j= max_j; j>F_max_j; j--) map_j2i[j-1] = FAR_FAR_AWAY;;
    i = F_max_i;
    j = F_max_j;
    *aln_score = F[i][j]; 
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
	    map_i2j [i-1] = FAR_FAR_AWAY;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = FAR_FAR_AWAY;
	    j--; 
	    break; 
	default: 
	    fprintf (stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}
/************************************************************/
/************************************************************/
int find_Calpha (Protein *protein, int  resctr, double ca[3]){

    Residue * res = protein->sequence + resctr;

    if ( ! res ) return 1;
    if ( ! res->Ca ) return 1;
    
    ca[0]= res->Ca->x;
    ca[1]= res->Ca->y;
    ca[2]= res->Ca->z;


    return 0;

}

/*************************************************************************/
/*************************************************************************/
int out_of_order_alignment ( Descr *descr1, Descr *descr2, Map *map, int *element_1_begin, int *element_1_end,
			     int *element_2_begin, int *element_2_end,
			     int longest_element_length, double ** similarity, double ** sim_in_element,
			     int *residue_map_i2j, int *residue_map_j2i, double * score_ptr) {


    int element_ctr_1, element_ctr_2;
    int resctr1, resctr2, ctr1, ctr2;
    int in_element_map_i2j[longest_element_length];
    int in_element_map_j2i[longest_element_length];
    double aln_score = 0;
    
    *score_ptr = 0.0;
    
    for (resctr1=0; resctr1< descr1->no_of_residues; resctr1++) {
	residue_map_i2j[resctr1] = FAR_FAR_AWAY;
    }

    for (resctr2=0; resctr2< descr2->no_of_residues; resctr2++) {
	residue_map_j2i[resctr2] = FAR_FAR_AWAY;
    }

    
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {
	
	element_ctr_2 = map->x2y[element_ctr_1];
	if (element_ctr_2 < 0) continue;

	ctr1 = 0; ctr2 = 0;
	for (resctr1=element_1_begin[element_ctr_1]; 
	     resctr1<= element_1_end[element_ctr_1]; resctr1++) {

	    ctr2 = 0;
	    for (resctr2=element_2_begin[element_ctr_2]; 
		 resctr2<= element_2_end[element_ctr_2]; resctr2++) {

		sim_in_element[ctr1][ctr2] = similarity[resctr1][resctr2];

		ctr2++;
	    }
	    ctr1++;  
	}

	smith_waterman_2 (ctr1, ctr2, sim_in_element, in_element_map_i2j,
			  in_element_map_j2i, &aln_score);
	*score_ptr += aln_score;

	
	ctr1 = 0; 
	for (resctr1=element_1_begin[element_ctr_1]; 
	     resctr1<= element_1_end[element_ctr_1]; resctr1++) {

	    residue_map_i2j[resctr1] = element_2_begin[element_ctr_2]+in_element_map_i2j[ctr1];
	    ctr1++;  
	}
	
	ctr2 = 0;
	for (resctr2=element_2_begin[element_ctr_2]; 
	     resctr2<= element_2_end[element_ctr_2]; resctr2++) {

	    residue_map_j2i[resctr2] = element_1_begin[element_ctr_1]+in_element_map_j2i[ctr2];
	    ctr2++;
	}
    }
    
    return 0;
}

/*************************************************************************/
/*************************************************************************/
int element_check (Protein *protein, Descr *descr) {

    int element_ctr;
    int * element_begin = protein->element_begin;
    int * element_end   = protein->element_end;

    
    printf (">>>>>>>>>>>>>>\n");
    printf (">>>>   %s \n", descr->name);
    for (element_ctr=0; element_ctr < descr->no_of_elements; element_ctr++) {
	printf (">>  %3d    %3d   %3d \n", element_ctr, element_begin[element_ctr], element_end[element_ctr]); 
    }
    printf (">>>>>>>>>>>>>>\n");

    return 0;

}




/*********************************************************************/
/*********************************************************************/
 
# if 0
	int residue_map_i2j_2[1000] = {0}, residue_map_j2i_2[1000] = {0};
	int i,j;
	
	for (j=0; j<no_res_2; j++) {
	    printf ("  %3d    %6d    -->   %6d \n", j, residue_map_j2i[j],  residue_map_j2i_2[j]);
	}
	exit(1);
# endif
