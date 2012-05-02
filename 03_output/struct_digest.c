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


int fill_scorecard (Representation *rep1, Representation *rep2, Map *map, Score *score) ;




int init_digest (Descr *qry_descr, Descr *tgt_descr, FILE ** digest_ptr) {


    FILE *digest = *digest_ptr;
    
    if (!digest) { /*open new one */
	char outname[MEDSTRING] = {'\0'};
	if (!options.outname[0] ) {
	    if ( qry_descr->name[0] &&  tgt_descr->name[0] ) {
		sprintf (outname, "%s_%s.struct_out", qry_descr->name, tgt_descr->name);
	    } else {
		sprintf (outname, "digest.struct_out");
	    }
	    sprintf (options.outname, "%s", outname);
	    
	} else {
	    sprintf (outname, "%s.struct_out", options.outname);
	}
	digest  = efopen (outname, "w");
    }
    
    if ( !digest) return 1;

    *digest_ptr = digest;
    
    if ( options.print_header) {
	fprintf ( digest,"%% columns: \n");
	fprintf ( digest,"%% query, target: structure names\n");
	fprintf ( digest,"%% geom_z:  z score for the orientational match \n");
	fprintf ( digest,"%% <dL>:    average length mismatch for matched SSEs \n");
	fprintf ( digest,"%% T:       total score assigned to matched SSEs \n");
	fprintf ( digest,"%% frac:    T divided by the number of matched SSEs \n");
	fprintf ( digest,"%% GC_rmsd: RMSD btw geometric centers "
		  "of matched SSEs (before postprocessing) \n");
	fprintf ( digest,"%% A:       (after postprocessing) the alignment score \n");
	fprintf ( digest,"%% aln_L:   (after postprocessing) the alignment length \n\n");
	fprintf ( digest,"%% %6s%6s %6s %6s  %6s %6s %6s %6s %6s %6s \n",
		  "query ", "target ", "geom_z", "<dL>", "  T  ", "frac",
		  "GC_rmsd", "rmsd  ", "A  ", "aln_L  " );	
    }

    return 0;

}


int write_digest(Descr *qry_descr, Descr *tgt_descr,
		 Representation * qry_rep, Representation * tgt_rep,
		 List_of_maps *list,FILE * digest) {


    if ( list ) {
	Score score;
	int rank_ctr, map_id;
	Map *current_map;
	for (rank_ctr=0; rank_ctr < list->no_maps_used &&
		 rank_ctr < options.number_maps_out; rank_ctr++) {
	    map_id = list->map_best[rank_ctr];
	    current_map = list->map+map_id;
	    fill_scorecard (tgt_rep, qry_rep, current_map, &score);
	    fprintf ( digest,
		      "%6s %6s %8.3lf %6.3lf %6.2lf %6.2lf %6.3lf %6.3lf %6.3lf %4d ",
		      qry_descr->name,
		      tgt_descr->name,
			      
		      score.z_score,
		      score.avg_length_mismatch,
			      
		      score.total_assigned_score,
		      score.fraction_assigned,
		  
		      score.rmsd,
			      
		      score.res_rmsd,
		      score.res_almt_score,
		      score.res_almt_length);
	    if (current_map->filename) fprintf (digest, "  %s", current_map->filename);
	    fprintf (digest, "\n");
	}
	
     } else {

	fprintf ( digest,
		  "%6s %6s %8.3lf %6.3lf %6.2lf %6.2lf %6.3lf %6.3lf %6.3lf %4d ",
		  qry_descr->name,
		  tgt_descr->name,
		  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0);

	fprintf (digest, "\n");
 
    }

    
    fflush  (digest);

    return 0;

}


int close_digest (clock_t CPU_time_begin, clock_t CPU_time_end, FILE *digest){

    fprintf (digest, "done   CPU:  %10.3lf s\n", (double)(CPU_time_end-CPU_time_begin)/CLOCKS_PER_SEC );
    fflush  (digest);
    fclose(digest);

    return 0;
}


/******************************************************************************************/
int fill_scorecard (Representation *rep1, Representation *rep2, Map *map, Score *score) {

    int sub_map_ctr;
    Map * current_map = NULL;
    int NX = rep1->N_full;
    int NY = rep2->N_full;
    
    if  (  !map)  return 1;
    
    memset (score, 0, sizeof(Score) );
    
    score -> z_score = map->z_score;
    score -> rmsd    = map->rmsd;
    score -> total_assigned_score = map->assigned_score;
	
    if ( score -> total_assigned_score >  0 ) {
	double avg_length_mismatch = 0;
	int i, j, map_size = 0;
	for (i=0; i < NX; i++ ) {
	    j =  map->x2y[i];
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

    score -> w_submap_z_score = map->z_score;
    sub_map_ctr = 0;
    while ( map->submatch_best && map->submatch_best[0] > -1 ) {

	current_map = map + current_map->submatch_best[0];
	if (  current_map->z_score >= 0 ) continue;
	score -> total_assigned_score += current_map->assigned_score;
	score -> w_submap_z_score     *=  current_map->z_score;
	score -> number_of_maps ++;
	if ( sub_map_ctr < 3) {
	    score->q_rmsd[sub_map_ctr] = quat_rmsd (map->q, current_map->q);
	}
	    
	sub_map_ctr++;
    }

    if ( NX > NY ) {
	score->fraction_assigned = score->total_assigned_score/NY;
    } else {
	score->fraction_assigned = score->total_assigned_score/NX;
    }
    score -> w_submap_z_score  = fabs(score -> w_submap_z_score);
	

    score->res_almt_score  = map->aln_score;
    score->res_almt_length = map->res_almt_length;
    score->res_rmsd        = map->res_rmsd;
	    	
    return 0;
}



