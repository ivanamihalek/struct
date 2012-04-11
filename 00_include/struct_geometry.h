#include "struct.h"

/* ss elements: */
# define HELIX 1
# define STRAND 2

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
    int type; /* helix or strand */
    char begin_id[PDB_HELIX_BEGIN_LEN+2]; /* this is a string identifier fomr PDB*/
    char end_id[PDB_HELIX_END_LEN+2]; 
    int begin, end; /* this may be added as post-processing step */
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
} Protein;


int protein_shutdown (Protein * protein);
