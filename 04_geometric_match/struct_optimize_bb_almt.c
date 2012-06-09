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

/**************************************************************/
/*  optimize the alignment on the backbone level, for each    */
/*  rotation+translation given in the map list                */
/**************************************************************/
int single_map_optimize_bb_almt (Protein * protein1, Protein * protein2, Map * map);

int optimize_backbone_alignment (Descr *descr1, Protein * protein1, Representation *rep1, 
				 Descr *descr2, Protein * protein2, Representation *rep2, 
				 List_of_maps *list){
    
     if ( list->no_maps_used == 0) return 1;
     
     int map_ctr,retval;
     
     Map *current_map;
       
     for (map_ctr=0; map_ctr<list->no_maps_used; map_ctr++) {
	 current_map = list->map+map_ctr;
	 retval = single_map_optimize_bb_almt (protein1, protein2, current_map);
	 if (retval) {
	     printf (" error optimize bb alignment   db:%s  query:%s \n", descr1->name, descr2->name);
	     exit (retval);
	 }
     }

     return 0;
}


int single_map_optimize_bb_almt (Protein * protein1, Protein * protein2, Map * map) {

    int step, no_steps = 10000;
    int reject, restart;
    double **R, T[3], T_new[3], T_best[3];
    double q[4], q_new[4], q_best[4], exp_s[4];
    double d0 = options.distance_tol_in_bb_almt;
    double best_score, new_score, score_diff;
    double dir_step_size = 0.01;
    double dice, total_score;
    double temperature = 1.0;
    double **similarity;
    
    if ( ! (R=dmatrix(3,3)) ) return 1;
    
    similarity = dmatrix (protein1->length, protein2->length);
    if ( !similarity ) return 1;
    
  
    /* the current best guess for the rotation and translation are in
       map->q (the rotation representeed as a 4-component quaternion; the components defined as doubles),
       and  map->T (3 component; double); to get the rotation matrix use 
       quat_to_R (q, R); defined in  04_geometric_match/struct_quaternion.c:30
       Map is defined in 00_include/struct.h:190      
    */
    
    printf ("\nin single_map_optimize_bb_almt():\n");
    printf ("initial alignment score: %8.4lf\n", map->aln_score);
    
    memcpy (&q[0], &(map->q[0]), 4*sizeof(double));
    memcpy (&T[0], &(map->T[0]), 3*sizeof(double));
    
    memcpy (&q_best[0], &(map->q[0]), 4*sizeof(double));
    memcpy (&T_best[0], &(map->T[0]), 3*sizeof(double));
    
    quat_to_R (q, R);
    
    closeness_score_for_bb_almt (map, protein1, protein2,
				 R, T, d0, similarity, &total_score);
    smith_waterman_2 (protein1->length, protein2->length, similarity,
		      map->x2y_residue_level, map->y2x_residue_level, &new_score);
    
    best_score = new_score;


    printf ("orig: %8.3lf\n", best_score);
    reject  = 0;
    restart = 0;
    for (step=1; step<=no_steps; step++) {

	memcpy (&(q_new[0]), &q[0], 4*sizeof(double));
	memcpy (&(T_new[0]), &T[0], 3*sizeof(double));
	
	/* pick a random move */
	dice = drand48();

	
# if 0
	if ( dice < 0.25) {
	    /* make  new  T*/
	    T_new[0] = T[0] + dir_step_size*(1-2*drand48());
	    T_new[1] = T[1];
	    T_new[2] = T[2];
	} else 	if ( dice < 0.5) {
	    T_new[0] = T[0];
	    T_new[1] = T[1] + dir_step_size*(1-2*drand48());
	    T_new[2] = T[2];	    
	} else 	if ( dice < 0.75) {
	    T_new[0] = T[0];
	    T_new[1] = T[1];
	    T_new[2] = T[2] + dir_step_size*(1-2*drand48());	    	   
	} else {
	    /* make  new q
	    random_q (exp_s, M_PI/360);
	    multiply (exp_s, q, 0, q_new);
	    quat_to_R (q_new, R);
	    memcpy (&(T_new[0]), &T[0], 3*sizeof(double));
	     */
	}
# endif	
	
	T_new[0] = T[0] + dir_step_size*(1-2*drand48());
	T_new[1] = T[1] + dir_step_size*(1-2*drand48());
	T_new[2] = T[2] + dir_step_size*(1-2*drand48());

	random_q (exp_s, M_PI/500);
	multiply (exp_s, q, 0, q_new);
	quat_to_R (q_new, R);
	
	/* redo Smith-Waterman */
	closeness_score_for_bb_almt (map, protein1, protein2,
				     R, T, d0, similarity, &total_score);
	smith_waterman_2 (protein1->length, protein2->length, similarity,
			  map->x2y_residue_level, map->y2x_residue_level, &new_score);
	
	/* MC step */	
	score_diff =  new_score -  best_score;
	if ( score_diff > 0   ||   drand48() < exp(score_diff/temperature)  ) {
	    memcpy (&(T[0]), &T_new[0], 3*sizeof(double));
	    memcpy (&(q[0]), &q_new[0], 4*sizeof(double));
	} else {
	    reject ++;
	}

	/* keep track of the best score we've found */
	if (  score_diff > 0 ) {
	    best_score = new_score;
	    memcpy (&(T_best[0]), &T_new[0], 3*sizeof(double));
	    memcpy (&(q_best[0]), &q_new[0], 4*sizeof(double));	    
	    printf ("step %4d    dice: %5.2lf     score: %8.4lf    diff: %8.4lf\n",
		    step, dice, new_score, score_diff);
	}


	if ( reject >= 100 ) {
	    memcpy (&(T_new[0]), &T_best[0], 3*sizeof(double));
	    memcpy (&(q_new[0]), &q_best[0], 4*sizeof(double));
	    restart++;
	}

	if (restart > 10) break;
	
    }
    
    /* after optimization replace map->q and map->T with the new values */
    memcpy (&(map->T[0]), &T_best[0], 3*sizeof(double));
    memcpy (&(map->q[0]), &q_best[0], 4*sizeof(double));

    map->aln_score =  best_score;
    
    
    free_dmatrix(R);
    free_dmatrix(similarity);

    return 0;
}


