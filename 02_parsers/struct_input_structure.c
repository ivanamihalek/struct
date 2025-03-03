/*
This source code is part of deconSTRUCT,
protein structure database search and backbone alignment application.
Written by Ivana Mihalek. Copyright (C) 2008-2025 Ivana Mihalek.

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

int pdb_input (FILE * fptr, char chain, Protein * protein,  Descr * description);
int db_input  (FILE * fptr, Descr * description);

int get_next_descr (int input_type, FILE * fptr,  char chain, Protein *protein, Descr * description) {

    int retval;
    
    if ( !fptr ) {
        fprintf (stderr, "Error in get_next_decr: read on closed fptr (?)\n");
        return -1;
    }

    if ( feof(fptr) ) return -1;

    
    if ( input_type == DB ) {
		retval = db_input (fptr, description);
    } else if (input_type == PDB ) {
		retval = pdb_input (fptr, chain, protein, description);
    } else {
		retval = -1;
    }
    
    return retval;

}

/********************************************************/
int pdb_input (FILE * fptr, char chain, Protein * protein, Descr * descr) {

    /* read in  the infor from the PDB file             */
    int retval = fill_protein_info (fptr, chain, protein);
    if (retval) return retval;

    /***********************************************************/
    if ( descr->db_file) { /* we read in the positions of SSEs */
		FILE *fptr  = NULL;
		int find_element_bounds (Protein *protein, Descr *descr) ;

		printf ("\nreading db from %s\n", descr->db_file);

		if ( !(fptr = efopen(descr->db_file, "r")) ) infox ("", 1);
		retval = db_input (fptr, descr);
		if ( retval ) return retval;

		find_element_bounds(protein, descr);
		protein->no_helices  = descr->no_of_helices;
		protein->no_strands  = descr->no_of_strands;
		protein->no_of_elements = descr->no_of_elements;
	
    /***********************************************************/
    }  else { /* we calculate the SSEs ourselves */

		/* find positions of SSEs on the sequence              */
		if (structure2sse (protein)) infox("Error finding SSEs.", 1);

		/* replace each SSE with a (directed) line:             */
		if (sse2descriptor (protein, descr))  infox("Error fitting lines to sse.",1);
    }

  
    return 0;
}


/*************************************************/
int db_input (FILE * fptr, Descr * descr) {

    int no_of_elements = 0;
    int no_of_helices  = 0, no_of_strands = 0;
    int retval;
    int element_ctr    = 0; 
    int descr_initialized =0;
    char line[BUFFLEN];
    SSElement * element;

    
    if (descr->element) descr_shutdown(descr);
    
    element_ctr = 0;
    memset (line,  0, BUFFLEN);
    
    while(fgets(line, BUFFLEN, fptr)!=NULL){

		if ( line[0] == '#' ) break; /* the end of record */
		if ( ! strncmp(line, "name", 4) ) {

		    sscanf (line, "%*s %s", descr->name);

		} else if ( ! strncmp(line, "number", 6) ) {

		    if ( strstr(line, "residues")) {
				sscanf ( line, "%*s %*s%*s  %d", &(descr->no_of_residues));
		    } else if ( strstr(line, "helices")){
				sscanf ( line, "%*s %*s %*s  %d", &no_of_helices);
				no_of_elements += no_of_helices;
				descr->no_of_helices   = no_of_helices;
				descr->no_of_elements += no_of_helices;
		    } else if (strstr(line, "strands")){
				sscanf ( line, "%*s %*s %*s  %d", &no_of_strands);
				/* if ( ! options.use_perp) no_of_strands/= 2; */
				no_of_elements += no_of_strands;
				descr->no_of_strands = no_of_strands;
				descr->no_of_elements += no_of_strands;
		    } else {
				fprintf (stderr, "Unrecognized record in the input:\n%s\n",
					 line);
				return 1;
		    }


		} else if ((!strncmp(line, "HELIX", 5)) || (!strncmp(line, "STRAND", 6))  ) {

		    if ( ! descr_initialized ) {
				descr_initialized = 1;
				if ( (retval = descr_init(descr) ) ) return retval;
		    }

		    /* store element info */
		    if ( element_ctr >= descr->no_of_elements ) {
			fprintf (stderr, "number of elements (%d) in %s\n",
				 element_ctr+1, descr->name) ;
			fprintf (stderr, "does not match the number declared (%d)\n",
				 descr->no_of_elements);
			return 1;
		    }
		    element = descr->element + element_ctr;

		    sscanf ( line, "%*s %d  %d  %s %s %lf %lf %lf %lf %lf %lf \n",
			     &(descr->element[element_ctr].type),
			     &(element->length), element->begin_id, element->end_id,
			     element->p, element->p+1, element->p+2,
			     element->cm, element->cm+1, element->cm+2);
		    descr->element[element_ctr].length = element->length;

		    element_ctr++;


		}
		memset (line,  0, BUFFLEN);
    }
    if ( !no_of_elements && feof(fptr) ) return -1;
    if ( descr->no_of_elements !=  element_ctr ) {
		fprintf (stderr, "number of elements (%d) in %s\n",
			 element_ctr, descr->name) ;
		fprintf (stderr, "does not match the number declared (%d)\n",
			 descr->no_of_elements);
		return 1;
    }
    
    /* make sure that p is normalized */    
    for ( element_ctr=descr->no_of_elements-1; element_ctr >=0; element_ctr-- ) {
		double norm, aux;
		int i;
		norm = 0;
		for (i=0; i<3; i++) {
		    aux = descr->element[element_ctr].p[i];
		    norm += aux*aux;
		}
		norm = sqrt (norm);
		for (i=0; i<3; i++) {
		    descr->element[element_ctr].p[i] /= norm;
		}
    }
    
    return 0;
}


/*************************************************/
int find_element_bounds (Protein *protein, Descr *descr) {
    
    int no_res = protein->length;
    int *element_begin, *element_end; /* "element" here means SSE */
    int *element_begin_pdb, *element_end_pdb;
    int resctr, element_ctr, num_pdb_id;
    char tmp[PDB_ATOM_RES_NO_LEN+1] = {'\0'};
    SSElement * element;
   
   
    if ( ! (element_begin     = emalloc (no_res*sizeof(int))) ) return 1;
    if ( ! (element_end       = emalloc (no_res*sizeof(int))) ) return 1;
    if ( ! (protein->sse_sequence = emalloc (protein->length*sizeof(int))) ) return 1;
    
    if ( ! (element_begin_pdb = emalloc (no_res*sizeof(int))) ) return 1;
    if ( ! (element_end_pdb   = emalloc (no_res*sizeof(int))) ) return 1;
  
    for (element_ctr=0; element_ctr < descr->no_of_elements; element_ctr++) {
		element = descr->element+element_ctr;
		/* get rid of the insertion tag, in case we
		haven't shaken it off along the way */
		memcpy (tmp,  element->begin_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
		element_begin_pdb[element_ctr] = atoi (tmp);
		
		memcpy (tmp,  element->end_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
		element_end_pdb[element_ctr] = atoi (tmp);
    }


    /* translate the beginning and the end of SSEs from pdb_tag to the position on the sequence*/
    element_ctr = 0;
    for (resctr=0; resctr<no_res; resctr++) {
	
	memcpy (tmp, protein->sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	num_pdb_id = atoi (tmp);
	
	/* break if we are outside of the last element */
	if ( num_pdb_id > element_end_pdb [descr->no_of_elements-1]) {
	    /*see for example k4j - the first residue has the number 399 */
	    continue;
	}
	
	if (  num_pdb_id == element_end_pdb[element_ctr] ) {
	    element_end [element_ctr] = resctr;
	    element_ctr++;
	}
	
	if ( num_pdb_id == element_begin_pdb[element_ctr] ) {
	    element_begin [element_ctr] = resctr;
	}
	/* also if there is some disagreement about where the first element starts: */
	if ( ! resctr && num_pdb_id >= element_begin_pdb[element_ctr] ) {
	    element_begin [element_ctr] = resctr;
	}
    }	      

    
    for (element_ctr=0; element_ctr < descr->no_of_elements; element_ctr++) {
		for (resctr = element_begin [element_ctr]; resctr <= element_end [element_ctr]; resctr++ ) {

			protein->sequence[resctr].belongs_to_helix  = 0;
			protein->sequence[resctr].belongs_to_strand = 0;

			if  ( descr->element[element_ctr].type == HELIX ) {
				protein->sequence[resctr].belongs_to_helix = element_ctr + 1;
				protein->sse_sequence[resctr] = HELIX;
			
			} else if   (descr->element[element_ctr].type == STRAND) {
				protein->sequence[resctr].belongs_to_strand = element_ctr + 1;
				protein->sse_sequence[resctr] = STRAND;
			}
		}
    }
    
    
    protein->element_begin     = element_begin;
    protein->element_end       = element_end;
 	
    free (element_begin_pdb);
    free (element_end_pdb);
    

    return 0;

}



/*************************************************/
int free_element (Descr * descr ) {

    SSElement * element = descr->element;

    free (element);
    descr-> no_of_elements = 0;

    
    return 0;
}


int alloc_element (Descr * descr) {

    int  no_of_elements = descr->no_of_elements;
    SSElement * element;
    
    element = emalloc ( no_of_elements*sizeof(SSElement) );
    if ( ! element) return 1;
    
    descr->element = element;
    descr->no_of_elements = no_of_elements;
    
    return 0;
}

/****************************************************************/
int descr_init ( Descr * descr ) {
		
    int retval;
    
    retval = alloc_element (descr);
    if (retval) return retval;
		
    return 0;
    
}


int descr_shutdown ( Descr * descr ) {
    
    int retval;
    
    if ( descr->element) {
	retval = free_element ( descr );
	if (retval) return retval;
    }

    memset (descr, 0, sizeof(Descr) );

    return 0;
}

/****************************************************************/
int check_input_type (FILE *fptr) {

    int input_type = 0;
    char line[LONGSTRING];
    
    if ( fgets(line,LONGSTRING,fptr)) {

	if ( ! strncmp (line, "name:",5) ) {
	    input_type = DB;
	} else if ( ! strncmp (line, "MODEL", 5) ) {
	    input_type = PDB;	    
	} else if ( ! strncmp (line, "TITLE", 5) ) {
	    input_type = PDB;	    
	} else if ( ! strncmp (line, "ATOM",  4) ) {
	    input_type = PDB;	    
	} else if ( ! strncmp (line, "HEADER", 6) ) {
	    input_type = PDB;	    
	} else if ( ! strncmp (line, "HETATM", 6) ) {
	    input_type = PDB;	    
	} else if ( ! strncmp (line, "REMARK", 6) ) {
	    input_type = PDB;	    
	}
    }
    rewind (fptr);

    return input_type;
}

