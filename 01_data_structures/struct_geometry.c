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

    int resctr, ectr;

    
    printf ("length:         %4d\n", protein->length);
    printf ("no_of_elements: %4d\n", protein->no_of_elements);
    printf ("no_helices:     %4d\n", protein->no_helices);
    printf ("no_strands:     %4d\n", protein->no_strands);

    for (ectr=0; ectr<protein->no_of_elements; ectr++) {
	printf ("\t %3d    %3d -- %3d",
		ectr, protein->element_begin[ectr],
		protein->element_end[ectr]);
	if (  protein->element_begin_pdb &&  protein->element_begin_pdb[ectr] ) {

	    printf ("\t  pdb: %3d -- %3d", protein->element_begin_pdb[ectr],
		    protein->element_end_pdb[ectr]);
	}
	
	printf ("\n");
   }

    printf ("\nSSE sequence\n");
    for (resctr=0; resctr<protein->length; resctr++) {
	printf ("%d", protein->sse_sequence[resctr]);
	if ( ! ((resctr+1)%50)) printf ("\n");
    }
    printf ("\n");

    printf ("\nhelices:  0x%x \n", protein->helix);
    if ( protein->helix ) {
	for (ectr=0; ectr<protein->no_helices; ectr++) {
	    SSElement * helix = protein->helix+ectr;
	    printf ("%3d    %4s  %4s  %c   %3d  %3d    %3d ", ectr,
		    helix->begin_id , helix->end_id, helix->chain,
		    helix->begin, helix->end, helix->length);
	    printf ("  %8.3lf %8.3lf %8.3lf ",
		    helix->p[0],  helix->p[1], helix->p[2]);
	    printf ("  %8.3lf %8.3lf %8.3lf ",
		    helix->cm[0],  helix->cm[1], helix->cm[2]);
	    printf ("\n");
	}
    }
    printf ("\nstrands:  0x%x \n", protein->strand);
    if ( protein->strand ) {
	for (ectr=0; ectr<protein->no_strands; ectr++) {
	    SSElement * strand = protein->strand+ectr;
	    printf ("%3d    %4s  %4s  %c   %3d  %3d    %3d ", ectr,
		    strand->begin_id , strand->end_id, strand->chain,
		    strand->begin, strand->end, strand->length);
	    printf ("  %8.3lf %8.3lf %8.3lf ",
		    strand->p[0],  strand->p[1], strand->p[2]);
	    printf ("  %8.3lf %8.3lf %8.3lf ",
		    strand->cm[0],  strand->cm[1], strand->cm[2]);
	    printf ("\n");
	}
    }

    
    printf ("\nAA  sequence\n");
    for (resctr=0; resctr<protein->length; resctr++) {
	Residue * residue = protein->sequence + resctr;
	printf ("%4d  %4s  %4s    %c  %c  h:%3d  s:%3d\n",
		resctr, residue->pdb_id, residue->res_type,
		residue->res_type_short, residue->chain,
		residue->belongs_to_helix, residue->belongs_to_strand);
		
    }
    printf ("\n");
   
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
