# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>
# include "nossa_pdb.h"
# include "nossa_utils.h"


# define  BUFFLEN 150

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25


typedef struct {
    char type [PDB_ATOM_ATOM_NAME_LEN+1];
    double x,y,z;
    int layer;
    int backbone;
} Atom;

# define  MAX_NO_ATOMS 40

typedef struct {
    char pdb_id[PDB_ATOM_RES_NO_LEN+1];
    char res_type[PDB_ATOM_RES_NAME_LEN+1];
    char res_type_short;
    int no_atoms;
    Atom  atom[MAX_NO_ATOMS];
    Atom * Ca;
} Residue;

typedef struct {
    int length;
    Residue * sequence;
} Protein;


/***********************************/


int almt_out (FILE *fptr, Protein * protein1, Protein *protein2,
	      int * residue_map_i2j, int * residue_map_j2i);
int read_pdb ( char * pdbname, Protein * protein, char chain);
int tfm_out  (FILE *fptr, Protein *protein1, Protein *protein2, int *residue_map_i2j);


# endif
