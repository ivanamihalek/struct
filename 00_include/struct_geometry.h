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

#include "struct.h"

/* ss elements: */
# define LOOP  1
# define HELIX 2
# define STRAND 4

# define MIN_CaS_IN_SSE 3  /* min number of Cas to be used for fitting of SSEs */


/* maximal allowed overlap between structural elements */
# define MAX_ALLOWED_OVERLAP 3

typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 100 /* so I can handle things like heme */

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+2];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    char chain;
    int no_atoms;
    Atom atom[MAX_NO_ATOMS];
    Atom *Ca;
    int interface;
    int solvent_accessible;
    int belongs_to_helix;
    int belongs_to_strand;
    int alt_belongs_to_helix; /* we allow a couple of residues overlap in SSEs */
    int alt_belongs_to_strand;
} Residue;



typedef struct {
    int type;                             /* helix or strand */
    char begin_id[PDB_HELIX_BEGIN_LEN+2]; /* this is a string identifier fomr PDB*/
    char end_id[PDB_HELIX_END_LEN+2];
    char chain;
    int begin, end;                       /* this may be added as post-processing step */
    int length;
    double p[3], cm[3];
} SSElement;


typedef struct {
    int length;
    Residue * sequence;
    int no_helices;
    SSElement *helix;
    int no_strands;
    SSElement *strand;
    int * sse_sequence;
    int * element_begin;
    int * element_end;
    int * element_begin_pdb;
    int * element_end_pdb;
    int no_of_elements;
} Protein;


int protein_shutdown (Protein * protein);
