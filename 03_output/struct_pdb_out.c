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

int pdb_output ( char *filename, double  **tfm_matrix, double * transl_vector, Residue * sequence, int no_res);

int transform_pdb (double  **tfm_matrix, double * transl_vector,
		   Residue * sequence, int no_res, Residue * sequence_new);


int write_tfmd_pdb ( Protein * tgt_protein, List_of_maps *list, Descr *tgt_descr, Descr *qry_descr) {
    
    if ( list->no_maps_used == 0) return 1;
     
    int rank_ctr, map_id, out_ctr;
    double ** R;
    char filename[LONGSTRING] = {'\0'};
    Map *current_map;
    Residue * sequence_new = emalloc (tgt_protein->length*sizeof(Residue));
    if (!sequence_new) return 1;
    if (!(R=dmatrix(3,3))) return 1; /* compiler is bugging me otherwise */


    out_ctr = 0;
	
    for (rank_ctr=0; rank_ctr<list->best_array_used && rank_ctr < options.number_maps_out; rank_ctr++) {
	
	map_id      = list->map_best[rank_ctr];
      	current_map = list->map+map_id;
	
	quat_to_R (current_map->q, R);

	/* copy the overall info about types, number of atoms etc from the original sequence*/ 
	memcpy (sequence_new, tgt_protein->sequence, tgt_protein->length*sizeof(Residue));
	/* transform the coordinates */ 
	transform_pdb (R, current_map->T, tgt_protein->sequence, tgt_protein->length, sequence_new);

	sprintf (filename, "%s.to_%s.%d.pdb", tgt_descr->name, qry_descr->name, out_ctr);
	if ( options.outdir[0])  {
	    char aux[MEDSTRING] = {'\0'};
	    sprintf (aux, "%s", filename);
	    memset (&filename[0], 0, LONGSTRING*sizeof(char));
	    sprintf (filename, "%s/%s",  options.outdir, aux);
	}
	pdb_output (filename, R, current_map->T, sequence_new, tgt_protein->length);
	out_ctr++;

	sprintf (current_map->filename, "%s", filename);
    }
    
    free_dmatrix(R);
    free(sequence_new);
    return 0;
}


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
int pdb_output ( char *filename, double  **tfm_matrix, double * transl_vector, Residue * sequence, int no_res){

    int i, j;
    int resctr, atomctr, serial;
    FILE * fptr;
    Atom * atomptr;
    
    /* open file */
    fptr = fopen ( filename, "w");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", filename);
	return 1;
    }
    /* output the trfm matrix, for the reference */
    fprintf (fptr, "REMARK \n" );
    fprintf (fptr, "REMARK  tfm matrix: \n" );
    fprintf (fptr, "REMARK \n" );
    for (i=0; i<3; i++) {
	fprintf (fptr, "REMARK ");
	for (j=0; j<3; j++) {
	    fprintf (fptr, "  %8.3lf ", tfm_matrix[i][j]);
	}
	fprintf (fptr, "  %8.3lf \n", transl_vector[i]);	
    }
    fprintf (fptr, "REMARK \n");
    fprintf (fptr, "REMARK \n");

    serial = 0;
    for ( resctr= 0; resctr < no_res; resctr ++ ) {
	for (atomctr=0; atomctr < sequence[resctr].no_atoms; atomctr++) {
	    serial ++;
	    atomptr = sequence[resctr].atom + atomctr;
	    fprintf (fptr,  "%-6s%5d  %-3s%1s%-3s %1c%4s%1s   %8.3lf%8.3lf%8.3lf\n", 
		     "ATOM",  serial ,  atomptr->type,   " ",   sequence[resctr].res_type,
		     sequence[resctr].chain, sequence[resctr].pdb_id,   " " ,
		     atomptr->x,  atomptr->y,   atomptr->z);
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
