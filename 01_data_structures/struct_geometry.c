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



/*************************************************************/
int protein_spit_out (Protein * protein) {

    printf ("length:         %4d\n", protein->length);
    printf ("no_of_elements: %4d\n", protein->no_of_elements);
    printf ("no_helices:     %4d\n", protein->no_helices);
    printf ("no_strands:     %4d\n", protein->no_strands);
    
    /*
    printf ("sequence: %d\n", );
    printf ("helix: %d\n", );
    printf ("strand: %d\n", );
    printf ("sse_sequence: %d\n", );
    printf ("element_begin: %d\n", );
    printf ("element_end: %d\n", );
    printf ("element_begin_pdb: %d\n", );
    printf ("element_end_pdb: %d\n", );
    */
   return 0;
}

/*************************************************************/
int protein_shutdown (Protein * protein) {

    if ( ! protein) return 0;
    
    protein->length = 0;
    if (protein->sequence) free(protein->sequence);
    
    if (protein->no_helices) protein->no_helices = 0;
    
    if (protein->helix) free (protein->helix);
    
    if (protein->no_strands) protein->no_strands = 0;
     
    if (protein->strand) free (protein->strand);
    if (protein->sse_sequence) free (protein->sse_sequence);

    
    return 0;
}
