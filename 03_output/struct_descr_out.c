# include "struct.h"

int descr_out (FILE * fptr, Descr * descr) {

    int element_ctr, i;
    SSElement * element;
    
    fprintf (fptr, "name: %s \n", descr->name);
    fprintf (fptr, "no of residues: %d \n", descr->no_of_residues);
    fprintf (fptr, "no of helices: %d \n", descr->no_of_helices);
    fprintf (fptr, "no of strands: %d \n", descr->no_of_strands);
    
    for (element_ctr=0; element_ctr < descr->no_of_elements; element_ctr++) {
	element = descr->element+element_ctr;
	if ( descr->element[element_ctr].type == HELIX ) {
	    fprintf (fptr, "HELIX  %3d %5d  %5s %5s ",
		     descr->element[element_ctr].type,
		     element->length,  element->begin_id, element->end_id);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->p[i]);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->cm[i]);
	    fprintf (fptr, "\n");
	} else {
	    fprintf (fptr, "STRAND %3d %5d  %5s %5s ",
		     descr->element[element_ctr].type,
		     element->length,  element->begin_id, element->end_id);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->p[i]);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->cm[i]);
	    fprintf (fptr, "\n");
	}
	
    }
    fprintf (fptr, "#\n"); /* end of record */
    
    return 0;
}


# if 0

/* junkyard */
	    
fprintf (fptr, "STRAND %3d %5d  %5s %5s ",
	 PERP,
	 element->length,  element->begin_id, element->end_id);
for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->perp[i]);
for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->foot[i]);
fprintf (fptr, "\n");

# endif
