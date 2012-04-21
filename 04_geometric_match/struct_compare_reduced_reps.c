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

/*****************************************************************************/
/*****************************************************************************/
int compare_reduced_reps (Representation *rep1, Representation *rep2,
		   List_of_maps *list, Score *score) {
    /* descr 1 is db or target, descr2 is query - important in postprocessing*/
    int retval;
    int map_ctr, best_ctr, i;
    int NX, NY, NX_eff, NY_eff;
    Map *current_map;

    /****************************************/
    /*  initialization                      */
    /****************************************/
    /*  shorthands   */
    NX     = rep1->N_full;
    NX_eff = rep1->N_compact;
    NY     = rep2->N_full;
    NY_eff = rep2->N_compact;
 
   
    if ( ! list->map )  return 1;
   
    
    /* the size has increased case */ 
    if ( NX > list->NX_allocated || NY > list->NY_allocated  ) {
	for ( map_ctr= 0; map_ctr< list->map_max; map_ctr++) {
	    if ( list->map) if ( free_map( list->map+map_ctr) ) return 1;
	    if ( initialize_map( list->map+map_ctr, NX, NY) ) return 1;
	}
	list->NX_allocated = NX;
	list->NY_allocated = NY;
    }

    for ( map_ctr= 0; map_ctr<list->map_max; map_ctr++) {
	(list->map+map_ctr)->x2y_size = NX;
	(list->map+map_ctr)->y2x_size = NY;
    }
    
    /****************************************/
    /* list  aliases                        */
    /****************************************/
    Map * map        = list->map;
    int map_max      = list->map_max;
    int best_max     = list->best_max;
    int *map_best    = list->map_best;


    /****************************************/
    /*  look for complementary maps         */
    /****************************************/
    map_ctr = 0;
  
    retval  = complement_match (rep1, rep2, map, map_max,
				&map_ctr, map_best, best_max, -1);
    if (retval) {
	return retval;
    }

    /****************************************/
    /*   output the best case stats         */
    /****************************************/
    memset (score, 0, sizeof(Score) );
    best_ctr = 0;
    while (  map_best[best_ctr] >  -1 ) {
	best_ctr ++;
    }
    
    if (best_ctr) {
	int sub_map_ctr;
	score -> z_score = map[map_best[0]].z_score;
	score -> rmsd = map[map_best[0]].rmsd;
	score -> total_assigned_score =  map[map_best[0]].assigned_score;
	
	if ( score -> total_assigned_score >  0 ) {
	    double avg_length_mismatch = 0;
	    int j, map_size = 0;
	    for (i=0; i < NX; i++ ) {
		j =  map[map_best[0]].x2y[i];
		if ( j <  0 )  continue;
		avg_length_mismatch += fabs(rep1->length[i] - rep2->length[j]);
		map_size++;
	    }
	    if (map_size) avg_length_mismatch /= map_size;
	    score -> avg_length_mismatch = avg_length_mismatch;
	    
	} else {
	    score -> avg_length_mismatch = 0;
	}
	score -> number_of_maps =1;

	current_map = map + map_best[0];
	score -> w_submap_z_score = map[map_best[0]].z_score;
	sub_map_ctr = 0;
	while ( current_map->submatch_best && current_map->submatch_best[0] > -1 ) {

	    current_map = map + current_map->submatch_best[0];
	    if (  current_map->z_score >= 0 ) continue;
	    score -> total_assigned_score += current_map->assigned_score;
	    score -> w_submap_z_score     *=  current_map->z_score;
	    score -> number_of_maps ++;
	    if ( sub_map_ctr < 3) {
		score->q_rmsd[sub_map_ctr] = quat_rmsd ( (map + map_best[0])->q, current_map->q);
	    }
	    
	    sub_map_ctr++;
	}

	if ( NX > NY ) {
	    score->fraction_assigned = score->total_assigned_score/NY;
	} else {
	    score->fraction_assigned = score->total_assigned_score/NX;
	}

 	score -> w_submap_z_score  = fabs(score -> w_submap_z_score);
	    
	
    } /* end if best_ctr */


    /****************************************/
    /*  aliases back to list                */
    /****************************************/
    list->map_max  = map_max;
    list->best_max = best_max;

   
    return 0;
    
}

