# include "struct.h"


int write_alignment (Protein *protein1, Protein *protein2,  List_of_maps * list) {

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
    char outname[BUFFLEN] = {'\0'};
    FILE *fptr;
    Residue *residue1, *residue2;
    
    int best_ctr = 0;
    int map_ctr  = list->map_best[best_ctr];
    Map *map    = list->map+map_ctr; /* we are still looking at the best map only */

    /* defined below: */
    int format_almt_line (Protein * protein, int * aligned_residues,
			 int first_pos, int last_pos, int chunk_size, char ret_str[]);
	
    
    if ( !(aligned_seq_x= emalloc ((protein1->length + protein2->length)*sizeof(int)))) return 1;
    if ( !(aligned_seq_y= emalloc ((protein1->length + protein2->length)*sizeof(int)))) return 1;
	
    last_mapped_x = map->x2y_residue_l_size-1;
    while (last_mapped_x>=0 && map->x2y_residue_level[last_mapped_x] < 0 ) last_mapped_x--;
    last_mapped_x++;

    last_mapped_y = map->y2x_residue_l_size-1;
    while (last_mapped_y>=0 && map->y2x_residue_level[last_mapped_y] < 0 ) last_mapped_y--;
    last_mapped_y++;
    
	
    first_mapped_x = -1;
    do {
	first_mapped_x ++;
    } while (first_mapped_x <map->x2y_residue_l_size
	     && map->x2y_residue_level[first_mapped_x] < 0);
	
    last_mapped_x = map->x2y_residue_l_size;
    do {
	last_mapped_x--;
    } while (last_mapped_x>=0 && map->x2y_residue_level[last_mapped_x] < 0 ); 



    first_mapped_y = -1;
    do {
	first_mapped_y ++;
    } while (first_mapped_y <map->y2x_residue_l_size
	     && map->y2x_residue_level[first_mapped_y] < 0);
	
    last_mapped_y = map->y2x_residue_l_size;
    do {
	last_mapped_y--;
    } while (last_mapped_y>=0 && map->y2x_residue_level[last_mapped_y] < 0 ); 



    i=first_mapped_x; j=first_mapped_y; pos = 0;
    while (i<map->x2y_residue_l_size ||  j<map->y2x_residue_l_size) {
	    
	if ( i == map->x2y_residue_l_size ) {
	    aligned_seq_x[pos] = -1;
	    aligned_seq_y[pos] =  j;
	    j++;
	    pos ++;
		
	} else if (j == map->y2x_residue_l_size) {
	    aligned_seq_x[pos] =  i;
	    aligned_seq_y[pos] = -1;
	    i++;
	    pos ++;
		
	} else if ( map->x2y_residue_level[i] >= 0  && map->y2x_residue_level[j] >= 0 ) {
	    aligned_seq_x[pos] =  i;
	    i++;
	    aligned_seq_y[pos] =  j;
	    j++;
	    pos ++;
		
	} else if ( map->x2y_residue_level[i] < 0 ) {
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

    /* as a temp measure: all maps have the  same name */
    sprintf (outname, "%s.alignment", options.outname);
    fptr  = efopen (outname, "w");
    if ( !fptr) exit (1);
    /* the plain old alignment */
    fprintf ( fptr, "%% alignment:\n%%\n");

    pos = 0;
    for ( chunk=0; chunk <= last_pos/chunk_size; chunk++) {

	format_almt_line (protein1, aligned_seq_x, pos, last_pos, chunk_size, formatted_almt_line);
	fprintf (fptr, "%% %s\n", formatted_almt_line);
	format_almt_line (protein2, aligned_seq_y, pos, last_pos, chunk_size, formatted_almt_line);
	fprintf (fptr, "%% %s\n", formatted_almt_line);
	fprintf (fptr, "%%\n%%\n");
	    
	pos += chunk_size;
    }

    /* find explicit rotation */
    if ( ! (R=dmatrix(3,3)) ) return 1; /* compiler is bugging me otherwise */
 
    quat_to_R (map->q, R);
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
	    
		point_rot_tr (ca1, R, map->T, rotated_ca1);
		
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
    free_dmatrix (R);

    free (aligned_seq_x);
    free (aligned_seq_y);

    fclose (fptr);

    

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
