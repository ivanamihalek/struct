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

int descr_out (FILE * fptr, Descr * descr) {

    int element_ctr, i;
    int my_fptr = 0;
    SSElement * element;

    if ( !fptr) {
	char filename[MEDSTRING] = {'\0'};
	sprintf (filename, "%s.db", descr->name);
	if ( ! (fptr=efopen(filename, "w")) )  return 1;
	my_fptr = 1;
    }
    
    fprintf (fptr, "name: %s \n", descr->name);
    fprintf (fptr, "number of residues: %d \n", descr->no_of_residues);
    fprintf (fptr, "number of helices: %d \n", descr->no_of_helices);
    fprintf (fptr, "number of strands: %d \n", descr->no_of_strands);
    
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

    if (my_fptr)  fclose(fptr);
	
    
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
