# include "struct.h"

int results_out (Descr *tgt_descr, Protein *tgt_structure, Representation * tgt_rep,
		 Descr *qry_descr, Protein *qry_structure, Representation * qry_rep,
		 List_of_maps *list, FILE * digest) {

    
    Score score;
    
    if (options.verbose) {
	write_maps (stdout, tgt_descr, qry_descr, list);
    }

    // score_map();
    if (options.postprocess) {
	write_tfmd_pdb  (tgt_structure, list, tgt_descr, qry_descr);
	write_alignment (tgt_structure, qry_structure, list);
    }
    
    write_digest(qry_descr, tgt_descr, digest, &score);               

    return 0;
}



/*********************************************************************************/
/*********************************************************************************/
# if 0
int score_map (Map *map, Score *score) {

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
	

	score->res_almt_score  = map->bb_aln_score;
	score->res_almt_length = map->res_almt_length;
	score->res_rmsd        = map->res_rmsd;


	
	    	
    } /* end if best_ctr */

    return 0;
}

# endif
