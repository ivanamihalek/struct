# include "struct.h"

int pdb_output ( char *filename, Residue * sequence, int no_res);
int transform_pdb (double  **tfm_matrix, double * transl_vector,
		   Residue * sequence, int no_res, Residue * sequence_new);


int write_tfmd_pdb ( Protein * tgt_protein, List_of_maps *list, Descr *tgt_descr, Descr *qry_descr) {
    
    /* for now,  we will just postprocess the  best map */
    int best_ctr = 0;
    int map_ctr = list->map_best[best_ctr];
    double ** R;
    char filename[MEDSTRING] = {'\0'};
    Map *map = list->map+map_ctr;
    Residue * sequence_new = emalloc (tgt_protein->length*sizeof(Residue));
    if (!sequence_new) return 1;
    
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
   
    quat_to_R (map->q, R);

    transform_pdb (R, map->T, tgt_protein->sequence, tgt_protein->length, sequence_new);

    sprintf (filename, "%s.rot_onto_%s.pdb", tgt_descr->name, qry_descr->name);
    pdb_output (filename, sequence_new, tgt_protein->length);
    free_dmatrix(R);
  
    return 0;
}


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
int pdb_output ( char *filename, Residue * sequence, int no_res){
    
    int resctr, atomctr, serial;
    FILE * fptr;
    Atom * atomptr;
    
    /* open file */
    fptr = fopen ( filename, "w");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", filename);
	return 1;
    }

    serial = 0;
    for ( resctr= 0; resctr < no_res; resctr ++ ) {
	for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	    serial ++;
	    atomptr = sequence[resctr].atom + atomctr;
	    fprintf (fptr,  "%-6s%5d  %-3s%1s%-3s %1s%4s%1s   %8.3lf%8.3lf%8.3lf\n", 
		"ATOM",  serial ,  atomptr->type,   " ",   sequence[resctr].res_type,
		"Z", sequence[resctr].pdb_id,   " " ,   atomptr->x,  atomptr->y,   atomptr->z);
	}
    }

    fclose (fptr);
    
    return 0;
}


/*******************************************************************************/
int transform_pdb (double  **tfm_matrix, double * transl_vector,
	       Residue * sequence, int no_res, Residue * sequence_new) {
    
    int resctr, atomctr;
    int i,j;
    double  x[3], new_x[3];
    
    for ( resctr= 0; resctr < no_res; resctr ++ ) {
	for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	    x[0] = sequence[resctr].atom[atomctr].x;
	    x[1] = sequence[resctr].atom[atomctr].y;
	    x[2] = sequence[resctr].atom[atomctr].z;
	    for (i=0; i <3; i++ ) {
		new_x[i] = 0.0;
		for (j=0; j <3; j++ ) {
		    new_x[i] += tfm_matrix[i][j]*x[j];
		}
		new_x[i] += transl_vector[i];
	    }
	    sequence_new[resctr].atom[atomctr].x = new_x[0];
	    sequence_new[resctr].atom[atomctr].y = new_x[1];
	    sequence_new[resctr].atom[atomctr].z = new_x[2];
	}
    }
    return 0;
}
