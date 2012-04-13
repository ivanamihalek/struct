# include "struct.h"

int pdb_input (FILE * fptr, char chain, Protein * protein,  Descr * description);
int db_input  (FILE * fptr, Descr * description);

int get_next_descr (int input_type, FILE * fptr,  char chain, Protein *protein, Descr * description) {

    int retval;
    
    if ( ! fptr ) {
	fprintf (stderr, "Error in get_next_decr: read on closed fptr (?)\n");
	return -1;
    }
    
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
    fill_protein_info (fptr, chain, protein);

    /* find positions of SSEs on the sequence           */
    if ( structure2sse (protein)) {
	fprintf ( stderr, "%s:%d: Error  finding SSEs.\n",
		  __FILE__, __LINE__);
	exit (1);
    }

    /* replace each SSE with a (directed) line:         */
    if ( sse2descriptor (protein, descr)) {
	fprintf ( stderr, "%s:%d: Error  fitting lines to sse.\n",
		  __FILE__, __LINE__ );
	exit (1);
    }

    /* fill in the remaining fields in the descriptor   */
    sprintf (descr->name, "%s", "anon");
    descr->no_of_residues = protein->length;

    return 0;
}


/*************************************************/
int db_input (FILE * fptr, Descr * descr) {

    int no_of_elements = 0;
    int no_of_helices  = 0, no_of_strands = 0;
    int type, retval;
    int element_ctr    = 0; 
    int descr_initialized =0;
    char line[BUFFLEN];
    SSElement * element;

    if ( feof(fptr) ) return -1;
    
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
		
	} else if (  (!strncmp(line, "HELIX", 5)) ||
		    (!strncmp(line, "STRAND", 6))  ) {
	    
	    if ( ! descr_initialized ) {
		descr_initialized = 1;
		if ( (retval = descr_init(descr) ) ) return retval;
	    }

	    
	    /* type check */
	    sscanf ( line, "%*s %d ", &type);
	    
	    /* if (  options.use_perp || */
            /* 		  (!options.use_perp  && type != PERP) ) { */
		
		
	    /* store element info */
	    if ( element_ctr >= descr->no_of_elements ) {
		fprintf (stderr, "number of elements (%d) in %s\n",
			 element_ctr+1, descr->name) ;
		fprintf (stderr, "does not match the number declared (%d)\n",
			 descr->no_of_elements);
		return 1;
	    }
	    element = descr->element + element_ctr;
	    /* the third field is an unused field (used to be sheet id, which was abandoned*/
	    sscanf ( line, "%*s %d  %*d   %d %s %s %lf %lf %lf %lf %lf %lf \n",
		     &(descr->element[element_ctr].type),
		     &(element->length), element->begin_id, element->end_id,
		     element->p, element->p+1, element->p+2, 
		     element->cm, element->cm+1, element->cm+2);
	    descr->element[element_ctr].length = element->length;

	    element_ctr++;
	    /*  } */
	    
	    
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
	} else if ( ! strncmp (line, "HEADER",6) ) {
	    input_type = PDB;	    
	} else if ( ! strncmp (line, "ATOM",4) ) {
	    input_type = PDB;	    
	} else if ( ! strncmp (line, "REMARK",6) ) {
	    input_type = PDB;	    
	}
    }
    rewind (fptr);

    return input_type;
}

