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


int write_alignment (Protein *protein1, Protein *protein2,  List_of_maps * list,
		     Descr *descr1,  Descr *descr2) {
				    

    int i,j;
    int pos, last_pos, chunk;
    int last_mapped_x, last_mapped_y;
    int first_mapped_x, first_mapped_y;
    int *aligned_seq_x, *aligned_seq_y;
    int chunk_size = 50;
    char formatted_almt_line[BUFFLEN]; /* BUFFLEN is 150 */
    double **R;
    int resctr_2, resctr_1;
    double ca1[3], ca2[3], rotated_ca1[3], d;
    char outname[LONGSTRING] = {'\0'};
    FILE *fptr;
    Residue *residue1, *residue2;
    
    int out_ctr  = 0;
    int rank_ctr = 0;
    int map_id = 0;
    Map * current_map;
 
    /* defined below: */
    int format_almt_line (Protein * protein, int * aligned_residues,
			  int first_pos, int last_pos, int chunk_size, char ret_str[]);
	
    
    if ( !(aligned_seq_x= emalloc ((protein1->length + protein2->length)*sizeof(int)))) return 1;
    if ( !(aligned_seq_y= emalloc ((protein1->length + protein2->length)*sizeof(int)))) return 1;

    if ( ! (R=dmatrix(3,3)) ) return 1; /* compiler is bugging me otherwise */
 



    out_ctr = -1;
	
    for (rank_ctr=0; rank_ctr<list->best_array_used && rank_ctr < options.number_maps_out; rank_ctr++) {
	
	map_id      = list->map_best[rank_ctr];
      	current_map = list->map+map_id;


	memset (aligned_seq_x, 0, (protein1->length + protein2->length)*sizeof(int));
	memset (aligned_seq_y, 0, (protein1->length + protein2->length)*sizeof(int));

 	last_mapped_x = current_map->x2y_residue_l_size-1;
	while (last_mapped_x>=0 && current_map->x2y_residue_level[last_mapped_x] < 0 ) last_mapped_x--;
	last_mapped_x++;

	last_mapped_y = current_map->y2x_residue_l_size-1;
	while (last_mapped_y>=0 && current_map->y2x_residue_level[last_mapped_y] < 0 ) last_mapped_y--;
	last_mapped_y++;
    
	
	first_mapped_x = -1;
	do {
	    first_mapped_x ++;
	} while (first_mapped_x < current_map->x2y_residue_l_size
		 && current_map->x2y_residue_level[first_mapped_x] < 0);
	
	last_mapped_x = current_map->x2y_residue_l_size;
	do {
	    last_mapped_x--;
	} while (last_mapped_x>=0 && current_map->x2y_residue_level[last_mapped_x] < 0 ); 



	first_mapped_y = -1;
	do {
	    first_mapped_y ++;
	} while (first_mapped_y <current_map->y2x_residue_l_size
		 && current_map->y2x_residue_level[first_mapped_y] < 0);
	
	last_mapped_y = current_map->y2x_residue_l_size;
	do {
	    last_mapped_y--;
	} while (last_mapped_y>=0 && current_map->y2x_residue_level[last_mapped_y] < 0 ); 



	i=first_mapped_x; j=first_mapped_y; pos = 0;
	while (i<current_map->x2y_residue_l_size ||  j<current_map->y2x_residue_l_size) {
	    
	    if ( i == current_map->x2y_residue_l_size ) {
		aligned_seq_x[pos] = -1;
		aligned_seq_y[pos] =  j;
		j++;
		pos ++;
		
	    } else if (j == current_map->y2x_residue_l_size) {
		aligned_seq_x[pos] =  i;
		aligned_seq_y[pos] = -1;
		i++;
		pos ++;
		
	    } else if ( current_map->x2y_residue_level[i] >= 0
			&& current_map->y2x_residue_level[j] >= 0 ) {
		aligned_seq_x[pos] =  i;
		i++;
		aligned_seq_y[pos] =  j;
		j++;
		pos ++;
		
	    } else if ( current_map->x2y_residue_level[i] < 0 ) {
		aligned_seq_x[pos] =  i;
		i++;
		aligned_seq_y[pos] = -1;
		pos ++;
		
		    
	    } else {
		aligned_seq_x[pos] = -1;
		aligned_seq_y[pos] =  j;
		j++;
		pos ++;
		
	    }
	    if ( i > last_mapped_x || j > last_mapped_y ) break;
	}

	last_pos = pos;

	/***************************************************************************/
	/* create the name for the output file and output it                       */
	out_ctr++;
	if (options.outdir[0] ) {
	    sprintf (outname, "%s/%s.%d.aln", options.outdir, options.outname,out_ctr );
	} else {
	    sprintf (outname, "%s.%d.aln", options.outname, out_ctr);
	}
	fptr  = efopen (outname, "w");
	if ( !fptr) exit (1);

	/***************************/
	/* the plain old alignment */
	fprintf ( fptr, "%% alignment  bewtween  %s  and  %s :\n%%\n", descr1->name, descr2->name);

	pos = 0;
	for ( chunk=0; chunk <= last_pos/chunk_size; chunk++) {

	    format_almt_line (protein1, aligned_seq_x, pos, last_pos, chunk_size, formatted_almt_line);
	    fprintf (fptr, "%% %s\n", formatted_almt_line);
	    format_almt_line (protein2, aligned_seq_y, pos, last_pos, chunk_size, formatted_almt_line);
	    fprintf (fptr, "%% %s\n", formatted_almt_line);
	    fprintf (fptr, "%%\n%%\n");
	    
	    pos += chunk_size;
	}

	/***************************/
	/* find explicit rotation  */
	quat_to_R (current_map->q, R);
	fprintf (fptr, "%% tfm matrix: \n" );
	for (i=0; i<3; i++) {
	    fprintf (fptr, "%%  ");
	    for (j=0; j<3; j++) {
		fprintf (fptr, "  %8.3lf ", R[i][j]);
	    }
	    fprintf (fptr, "  %8.3lf \n", current_map->T[i]);
	}

	
	/* mapped secondary structure */
	fprintf (fptr, "%%\n");
	fprintf (fptr, "%% mapped elements of secondary structure \n" );
	fprintf (fptr, "%% SSE1-->SSE2   type   cosine  exp_score   chain1  range1  chain2  range2 \n"); 
	
	print_map (fptr, current_map, descr1, descr2,   protein1, protein2, 0);
	fprintf (fptr, "%%\n");
	fprintf (fptr, "%%\n");
	
	

	/* mapped Ca's structure */
	fprintf (fptr, "%%\n");
	fprintf (fptr, "%% distances of the aligned Ca's:\n%%\n" );
	fprintf (fptr, "%% pdb_id_1  type_1    distance  pdb_id_2  type_2\n"); 
	for (pos=0; pos < last_pos; pos++) {
	    resctr_1 = aligned_seq_x[pos];
	    resctr_2 = aligned_seq_y[pos];
	
	    residue1 = protein1->sequence+resctr_1;
	    residue2 = protein2->sequence+resctr_2;

	    if ( resctr_1 < 0) {
		fprintf (fptr, "%7c  %7c ",  '-', '-');
	    } else {
		fprintf (fptr, "%7s  %7c ", residue1->pdb_id, residue1->res_type_short);
	    }
	    if ( resctr_1 < 0 || resctr_2 < 0 ) {
		fprintf (fptr, " %10c ", '-');
	    } else {
		/** evaluate the distance */

		if (  residue1->Ca && residue2->Ca) {
	    
		    ca1[0] = residue1->Ca->x;
		    ca1[1] = residue1->Ca->y;
		    ca1[2] = residue1->Ca->z;
	    
		    ca2[0] = residue2->Ca->x;
		    ca2[1] = residue2->Ca->y;
		    ca2[2] = residue2->Ca->z;
	    
		    point_rot_tr (ca1, R, current_map->T, rotated_ca1);
		
		    d = two_point_distance (rotated_ca1, ca2);
		    fprintf (fptr, " %10.2lf ", d);

		} else { /* survive cases when Ca not given */
		    fprintf (fptr, " NA ");
		}
		     
	    }
	    if ( resctr_2 < 0) {
		fprintf (fptr, "%7c  %7c ",  '-', '-');
	    } else {
		fprintf (fptr, "%7s  %7c ", residue2->pdb_id,  residue2->res_type_short);
	    }
	    fprintf (fptr, "\n");
	}
	fclose (fptr);

    }

    
    free_dmatrix (R);
    
    free (aligned_seq_x);
    free (aligned_seq_y);


    return 0;
}
/*********************************************/

int format_almt_line (Protein * protein, int * aligned_residues,
		      int first_pos, int last_pos, int chunk_size, char ret_str[]) {

    int pos, i, resctr;

    memset (ret_str, 0 , BUFFLEN*sizeof(char));

    
    /* look for the first label in this line*/
    pos = first_pos;
    resctr = -1;
    for (i=0; i< chunk_size && pos < last_pos; i++) {
	resctr = aligned_residues[pos];
	if ( resctr > -1) break;
	pos ++;
    }
    if (  resctr > -1 ) {
	sprintf (ret_str, " %3s ", protein->sequence[resctr].pdb_id);
    } else {
	sprintf (ret_str, " %3s ", "");
    }
    
    pos = first_pos;	    
    /* the next  chunk_size cahracters in the  protein */
    for (i=0; i< chunk_size && pos < last_pos; i++) {
	resctr = aligned_residues[pos];
	if (  resctr > -1 ) {
	    sprintf (ret_str,  "%s%c", ret_str, protein->sequence[resctr].res_type_short);
	} else {
	    sprintf (ret_str,  "%s%c", ret_str, '.');
	}
	pos++;

    }


    /* look for the last label*/
    pos -- ;
    resctr = -1;
    for (i=0; i< chunk_size && pos > 0; i++) {
	resctr = aligned_residues[pos];
	if ( resctr > -1) break;
	pos --;
    }
    if (  resctr > -1 ) {
	sprintf (ret_str, "%s %3s ", ret_str, protein->sequence[resctr].pdb_id);
    }
	    
 
    return 0;

}
