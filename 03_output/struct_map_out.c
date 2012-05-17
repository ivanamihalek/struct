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

/************************************************************************************/
/************************************************************************************/
/************************************************************************************/

/* output for a human reader (in verbose option)*/

int write_maps (FILE * fptr, Descr *descr1, Descr *descr2, List_of_maps *list) {

    int map_id, rank_ctr;
    int recursive_map_out (Map * map,  Descr * descr1, Descr * descr2, 
			   Protein *protein1, Protein *protein2, int depth);
    
    for (rank_ctr=0; rank_ctr<list->no_maps_used && rank_ctr < options.number_maps_out; rank_ctr++) {
	map_id = list->map_best[rank_ctr];
	recursive_map_out (list->map+map_id, descr1, descr2, NULL, NULL, 0);
    }

    return 0;
}

/************************************************************************************/
int recursive_map_out (Map * map, 
		       Descr * descr1, Descr * descr2,
		       Protein *protein1, Protein *protein2,
		       int depth) {
    
    int  ctr;
    
    if ( ! depth) {
	printf ("\n####################################\n");
	printf ("**  %s   %s ", descr1->name, descr2->name);
	if (map->submatch_best && map->submatch_best[0] > -1) {
	    printf ("  with submap (%d)", map->submatch_best[0]+1);
	}
	printf ("\n");
    } else {
	for (ctr=0; ctr < depth; ctr++) printf ("\t");
    }
    
    printf ("   geometric z_score: %6.3lf ",  map->z_score);
    if ( options.postprocess )  {
	 printf ("   alignment score:  %6.3lf ", map->aln_score);
    }
    printf ("\n");
    
    print_map (stdout, map, descr1, descr2,   protein1, protein2,  depth);
    
    if (map->submatch_best && map->submatch_best[0] > -1 ) {
	
	int sub_best_ctr, sub_map_ctr;
	sub_best_ctr = 0;
	while ( sub_best_ctr < options.number_maps_cpl &&
		(sub_map_ctr = map->submatch_best[sub_best_ctr] ) > -1 ) {
	    for (ctr=0; ctr < depth; ctr++) printf ("\t");
	    printf ("map %3d   submatch  complementarity z-score:  %6.3lf \n",
		    sub_map_ctr+1, map->compl_z_score); 
	    recursive_map_out (map+sub_map_ctr, descr1, descr2,
			       protein1, protein2, depth+1); 
	    sub_best_ctr ++;
	    if (sub_best_ctr==1) break; /*output only the best */
	} 
    }
    return 0; 
}

/************************************************/
/************************************************/
/************************************************/
int print_transformation (FILE* fptr, Map * map, int depth){

     
    double **R;
    int i, j, ctr;
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    quat_to_R (map->q, R);
    //fprintf (fptr,"transformation\n");
    for (i=0; i<3; i++) {
	if ( depth) for ( ctr=0; ctr < depth; ctr++) fprintf (fptr, "\t");
	for (j=0; j<3; j++) fprintf (fptr,"%8.3lf ", R[i][j]);
	fprintf (fptr,"%8.3lf  \n", map->T[i]);
    }

    free_dmatrix (R);
   
  
    return 0;
    
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
int print_transformation_0 (FILE* fptr, Map * map, int depth){

    double **R;
    int i, j, ctr;
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    quat_to_R (map->q, R);
    fprintf (fptr,"rotation\n");
    for (i=0; i<3; i++) {
	if ( depth) for ( ctr=0; ctr < depth; ctr++) fprintf (fptr, "\t");
	for (j=0; j<3; j++) fprintf (fptr,"%8.3lf ", R[i][j]);
	fprintf (fptr,"\n");
    }

    free_dmatrix (R);
   
    return 0;
    
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
int print_map_for_postproc (FILE *fptr, Map * map,
			    Representation *X_rep,  Representation *Y_rep) {
    
    int i, j, k;

    fprintf ( fptr, "translations\n");
    for (i=0; i<map->x2y_size; i++) {
	    
	j = map->x2y[i];
	if ( j < 0 ) continue;
	fprintf ( fptr, " %3d ", i+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  X_rep->translation[i][k]);
	}
	fprintf ( fptr, " %3d ", j+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  Y_rep->translation[j][k]);
	}
	fprintf ( fptr, "\n");
    }

    fprintf ( fptr, "cms\n");
    for (i=0; i<map->x2y_size; i++) {
	    
	j = map->x2y[i];
	if ( j < 0 ) continue;
	fprintf ( fptr, " %3d ", i+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  X_rep->cm[i][k]);
	}
	fprintf ( fptr, " %3d ", j+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  Y_rep->cm[j][k]);
	}
	fprintf ( fptr, "\n");
    }
    fprintf (fptr,"origin1 ");
    for (k=0; k<3; k++) {
	fprintf (fptr,"%8.3lf ", X_rep->origin[k]);
    }
    fprintf ( fptr, "\n");
    fprintf (fptr,"origin2 ");
    for (k=0; k<3; k++) {
	fprintf (fptr,"%8.3lf ", Y_rep->origin[k]);
    }
    fprintf ( fptr, "\n");
    
    return 0;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
int print_map (FILE *fptr, Map * map, Descr * descr1, Descr * descr2,
	       Protein * protein1, Protein *protein2, int tab) {
    
    int i, index_x, index_y;
    char type;
    
    /* header */
    if (options.print_header) {
	fprintf ( fptr, "%%");
	if (tab) fprintf ( fptr, "\t");
	fprintf ( fptr, " %9s   %4s %6s   %8s  %9s  %9s\n", "SSE map", "type",
		  " cosine", "Gauss wt.", "res range", "res range");
	if ( descr1 && descr2) {
	    fprintf ( fptr, "%%");
	    if (tab) fprintf ( fptr, "\t");
	    fprintf ( fptr, "%5s-->%5s   %4s   %6s %8s in %5s   in %5s\n",
		      descr1->name, descr2->name,
		      "", "", "",  descr1->name, descr2->name);
	}
    }
    /* the actual map */
    for (i=0; i<map->x2y_size; i++) {
	    
	index_x =  i;
	index_y =  map->x2y[i];
	if ( index_y <  0 ) continue;

	type = ' ';
	if ( descr1 ) {
	    if (descr1->element[index_x].type == HELIX) {
		type = 'H';
	    } else {
		type = 'S';
	    }
	}
	
	fprintf ( fptr, "%%");
	if (tab) fprintf ( fptr, "\t");
	fprintf ( fptr, " %2d --> %2d  %4c  %6.2lf   %8.1le ",
		  index_x+1, index_y+1,  type, map->cosine[index_x][index_y],
		map->image[index_x][index_y]);
	if ( descr1 && descr2) {
	    char c1, c2;
	    if ( descr1->element[index_x].chain == '\0' ||
		 descr1->element[index_x].chain == ' ') {
		c1 = '-';
	    } else {
		c1 =  descr1->element[index_x].chain;
	    }
	    if ( descr2->element[index_y].chain == '\0' ||
		 descr2->element[index_y].chain == ' ') {
		c2 = '-';
	    } else {
		c2 =  descr2->element[index_y].chain;
	    }

	    
	
	    fprintf (fptr,  "    %c %4s %4s   %c %4s %4s",
		     c1, descr1->element[index_x].begin_id,
		     descr1->element[index_x].end_id,
		     c2, descr2->element[index_y].begin_id,
		     descr2->element[index_y].end_id);
	}
	fprintf ( fptr, "\n");
    }
    fprintf ( fptr, "%%\n");
    fprintf ( fptr, "%%\n");
    fprintf ( fptr, "%%");
    if (tab) fprintf ( fptr, "\t");
    fprintf ( fptr, "total assigned score %8.4lf \n", map->assigned_score);

     
    return 0;
}



