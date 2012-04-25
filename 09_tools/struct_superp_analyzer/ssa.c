# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include "pdb.h"
# include "utils.h"

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define  LONGSTRING  250
# define  MEDSTRING  100
# define  SHORTSTRING  25
# define ALMT_NAME_LENGTH 30
# define BUFFLEN  250
/* pdb parser errors */
# define ERR_SSE_NONE    2
# define ERR_CA_TRACE    4
# define ERR_NONSENSE    8
# define ERR_CONTAINER  16
# define ERR_OVERLAP    32
# define ERR_MAX_ATOMS  64
# define ERR_NO_FILE_OR_CHAIN  128

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int * seq_gaps;
    int * column_gaps;
}  Alignment;

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
    int length;
    Residue * sequence;
    int no_helices;
    int no_strands;
    int * sse_sequence;
    int *prot2almt;
    int *almt2prot;
} Protein;

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

char single_letter ( char code[]);
int  read_pdb ( char * pdbname,  char chain, Protein * protein);
int  string_clean ( char* string, int length);

int tokenize ( char token[MAX_TOK][MEDSTRING], int * max_token, char * line , char comment_char);
FILE * efopen(char * name, char * mode);
int read_clustalw ( char * cwname, Alignment * alignment);
int  struct_almt_mapping (Protein * protein, Alignment * alignment, int sequence_id);
double  two_point_distance (double point1[3], double point2[3]);

int main ( int argc, char * argv[]) {

    int s, retval, almt_pos;
    Alignment alignment;
    Protein * protein;
    
    
    if ( argc <2 ) {
	fprintf ( stderr, "Usage: %s  <msf file>\n", argv[1]);
	exit (1);
    }

    /* read in the alignment */
    read_clustalw ( argv[1], &alignment);

    /* assign space for as amny structures as there are names in the alignment */
    protein = emalloc ( alignment.number_of_seqs*sizeof(Protein) );
    /*  for each name in the alignment read in the pdb fileW */
    for ( s=0; s<alignment.number_of_seqs; s++ ) {
	//printf ("%s\n", alignment.name[s] );

	/* read in the pdb file */
	if ( read_pdb ( alignment.name[s],  '\0', protein+s) ) return 1;
    
	/* map the pdb to the alignment */
	if ( ! (protein[s].prot2almt = (int *) emalloc (protein[s].length*sizeof(int))) ) exit (1);
	if ( ! (protein[s].almt2prot = (int *) emalloc (alignment.length*sizeof(int))) )exit (1);
	retval    = struct_almt_mapping (protein+s, &alignment, s);
	if (retval) {
	    fprintf (stderr, "Error mapping %s.\n", alignment.name[s]);
	    exit(retval);
	}
    }
    
    /*for each position in the alignment the alignment score */
    int s1, s2, no_points, pos;
    double score   = 0.0;
    double score_s = 0.0;
    double ca1[3], ca2[3];
    double dist, d0 = 10.0;
    
    for (almt_pos=0; almt_pos < alignment.length; almt_pos++ ) {
	no_points = 0;
	score_s = 0;
	
	for ( s1=0; s1<alignment.number_of_seqs; s1++ ) {
	    pos = protein[s1].almt2prot[almt_pos];
	    if ( ! protein[s1].sequence[pos].Ca) continue;
	    ca1[0] = protein[s1].sequence[pos].Ca->x;
	    ca1[1] = protein[s1].sequence[pos].Ca->y;
	    ca1[2] = protein[s1].sequence[pos].Ca->z;
	    
	    for ( s2=s1+1; s2<alignment.number_of_seqs; s2++ ) {
		pos = protein[s2].almt2prot[almt_pos];
		if ( ! protein[s2].sequence[pos].Ca) continue;
		ca2[0] = protein[s2].sequence[pos].Ca->x;
		ca2[1] = protein[s2].sequence[pos].Ca->y;
		ca2[2] = protein[s2].sequence[pos].Ca->z;
		
		dist =two_point_distance(ca1, ca2);
		score_s += exp(-dist/d0);
		no_points++;
	    }
	}
	if (no_points) score += score_s/no_points;
    }

    printf (" score:  %8.3lf  %8.3lf   \n", score, score/alignment.length);
    return 0;
}

/*************************************************************************/
double  two_point_distance (double point1[3], double point2[3] ) {
    int i;
    double aux;
    double d = 0;
    for (i=0; i<3; i++) {
	aux = point1[i] - point2[i];
	d  += aux*aux;
    }
    return  sqrt (d);
}

/************************************************************/
int count_gaps (Alignment * alignment) {

    int s, c;
    alignment->seq_gaps    = (int *) emalloc (alignment->number_of_seqs*sizeof(int));
    if (!alignment->seq_gaps) return 1;
    alignment->column_gaps = (int *) emalloc (alignment->length*sizeof(int));
    if (!alignment->column_gaps) return 1;
    for ( s=0; s<alignment->number_of_seqs; s++ ) {
	for ( c=0; c<alignment->length; c++) {
	    if ( alignment->sequence[s][c] == '.' ) {
		alignment->column_gaps[c] ++;
		alignment->seq_gaps[s] ++;
	    }
	}
    }
    return 0;
}
/************************************************************/

/*****************************************/
int read_clustalw ( char * cwname, Alignment * alignment){
    
    FILE * fptr = NULL;
    char line[BUFFLEN];
    int  number_of_seqs, almt_length, ctr;
    int * seq_pos, pos;
    char ** sequence;
    char ** name;
    char curr_name[BUFFLEN];
     
    /* open file */
    fptr = efopen ( cwname, "r");
    if ( !fptr ) return 1;
    
    /* find the alignment length info */
    almt_length = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( strstr(line, "MSF:" ) ){
	    sscanf (line, "%*s %d", &almt_length);
	    break;
	}
    }
    if ( almt_length ) {
	/* printf ( "Alignment length in %s is %d.\n", cwname, almt_length); */
    } else {
	fprintf ( stderr, "Alignment length info not found in %s. Is the format gcg?\n", cwname);
	return 1;
    }

    /* determine the number of sequences */
    number_of_seqs = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( ! strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) number_of_seqs++;
    }
    if ( number_of_seqs ) {
	/* printf ( "Number of sequences in %s is %d.\n", cwname, number_of_seqs); */
    } else {
	fprintf ( stderr, "No sequences found in %s. Is the format gcg?\n", cwname);
	return 1;
    } 
    
    /* allocate */
    sequence = chmatrix (number_of_seqs, almt_length);
    if ( !sequence ) return 1;
    name     = chmatrix (number_of_seqs, ALMT_NAME_LENGTH);
    if ( !name ) return 1;
    seq_pos = (int *) calloc ( number_of_seqs, sizeof(int));
    if ( !seq_pos ) return 1;
    
    /* read in */
    rewind(fptr);
    ctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if (!  strncmp (line, "//", 2) ) break;
	if ( strstr(line, "Name:" ) ) {
	    sscanf (line, "%*s %s", name[ctr]);
	    ctr ++;
	}
    }
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( isspace (line[0] ) ) continue;
	sscanf (line, "%s", curr_name);
	ctr = 0;
	while (  ctr <number_of_seqs &&  strcmp (name[ctr], curr_name) ) ctr++;
	if ( ctr >= number_of_seqs ) {
	    fprintf ( stderr, "The name %s not found in the header of %s.\n", curr_name, cwname);
	    return 1;
	}
	pos = 0;
	while ( ! isspace(line[pos]) ) pos++;
	while  (line[pos] != '\n' && pos < BUFFLEN) {
	    if ( !  isspace(line[pos] ) ){
		if ((line[pos]>=97)&&(line[pos]<=122)) {line[pos] -= 32;} /* --> turn to uppercase */
		if ( line[pos]==126)                   {line[pos]  = 46;} /* turn tweedle to dot */
		sequence [ctr] [ seq_pos[ctr] ] = line[pos];
		seq_pos[ctr]++;
	    }
	    pos ++;
	}
    }
    fclose(fptr);

    /* sanity check */
    for (ctr=0; ctr < number_of_seqs; ctr++ ) {
	if ( seq_pos[ctr] != almt_length ) {
	    fprintf (stderr, "Sequence %s is shorter (%d position) than the alignment.\n", name[ctr],  seq_pos[ctr]);
	    return 1;
	}
    }

    /* return values */
    alignment->number_of_seqs = number_of_seqs;
    alignment->length         = almt_length;
    alignment->sequence       = sequence;
    alignment->name           = name;

    /* free */
    free (seq_pos);

    { 
	int count_gaps (Alignment * alignment);
	count_gaps (alignment);
    }
    return 0;
}


/*****************************************************************/
int read_pdb ( char * pdbname,  char chain, Protein * protein) {
    
    int retval;
    FILE * fptr = NULL;
    int fill_protein_info ( FILE * fptr,  char chain, Protein * protein);
    char filename[MEDSTRING];
    
    sprintf (filename, "%s.pdb", pdbname);
    /* open file */
    fptr = fopen ( filename, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", filename);
	return ERR_NO_FILE_OR_CHAIN;
    }

    retval = fill_protein_info (fptr, chain, protein);
    if (retval) return retval;
    
    fclose(fptr);
    return 0;
}

/*****************************************************************/
int fill_protein_info ( FILE * fptr,  char chain, Protein * protein) {

    /* TODO for the moment we rely on PDB annotation
       to extract structural elements - that should be changed */
    Residue * sequence = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+2];
    /* res name: 4 digits + insertion code + \0 */
    char oldrestype [PDB_ATOM_RES_NAME_LEN+2];
    char tmp[BUFFLEN], *auxptr;
    char atomtypes_read_in[BUFFLEN];
    char old_chain;
    int atomctr, resctr,  no_res,ctr, nonblank;
    int retval;
    int chain_found;
    int ca_trace;
    char single_letter ( char code[]);
    
    int has_backbone (Residue * sequence, int from, int to);
    
    /********************************************/
    /********************************************/
    /* cleanup                                  */
    memset (protein, 0, sizeof(Protein) );

    /********************************************/
    /********************************************/
    /* count residues                           */
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    resctr = 0;
    chain_found = 0;
    old_chain = '\0';
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	if (resctr) {
	    if ( ! strncmp(line,"END", 3) ||  (chain && line[PDB_ATOM_CHAINID] != old_chain) )
		break;
	}
	if (chain  && line[PDB_ATOM_CHAINID] != chain) continue;
	chain_found  = 1;
	
	if( ! strncmp(line,"ATOM", 4)){
	    
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
		/* handling the case when the chain is not given, meaning: "take the first chain" */ 
		old_chain = line[PDB_ATOM_CHAINID];
		resctr ++;
	    }
	} 
    }

    /* sanity: */
    if ( chain && ! chain_found) {
	fprintf (stderr, "Chain %c not found.\n", chain);
	return ERR_NO_FILE_OR_CHAIN;
    }

    no_res = resctr;
    
    if ( !no_res ) return -1;  /* take it as the end of the read */
    
    /* allocate space */
    sequence = emalloc ( no_res*sizeof (Residue));
    if ( ! sequence ) return 1;
  
    
    /*********************************************/
    /*********************************************/
    /*   read in residue numbers and atom coords */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    old_chain = '\0';

    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    resctr= -1;
    atomctr = 0;
    ca_trace = 1;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	
	if ( resctr > -1) {
	    if  (! strncmp(line,"END", 3)  ||  (chain && line[PDB_ATOM_CHAINID] != old_chain))
	    break;
	}
	if ( chain  && line[PDB_ATOM_CHAINID] != chain ) continue;
	
	if( ! strncmp(line,"ATOM", 4) ){
 	   /* if it's a hydrogen - skip */
	    if ( line[PDB_ATOM_ATOM_NAME] == 'H'
		 ||  line[PDB_ATOM_ATOM_NAME+1] == 'H') continue;
	    /* adjust the counters */ 
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		/*+1 in  PDB_ATOM_RES_NO_LEN+1 means I am including the insertion code
		  in the identifier */
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		strncpy (oldrestype, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN+1]   = '\0';
		oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
		
		/* handling the case when the chain is not given, meaning: "take the first chain" */ 
		old_chain = line[PDB_ATOM_CHAINID];
		
		resctr ++;
		if ( resctr >= no_res ) {
		    fprintf (stderr, "Error reading pdb: resctr:%d   no res: %d\n",
			     resctr, no_res);
		    return ERR_NONSENSE;
		}
		atomctr = 0;
		/* keep track of atom types we have read in */
		memset (atomtypes_read_in, 0, BUFFLEN*sizeof(char));
		atomtypes_read_in[0] = '_';
		
		sequence[resctr].no_atoms = 1;
		strncpy ( sequence[resctr].pdb_id, oldresno, PDB_ATOM_RES_NO_LEN+2);
		sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN+1]   = '\0';
		
		strncpy ( sequence[resctr].res_type, oldrestype, PDB_ATOM_RES_NAME_LEN+1);
		sequence[resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';
		sequence[resctr].res_type_short  =
		    single_letter ( sequence[resctr].res_type );
		/* modified residues are ok for the purposes here */
		/* unless they are sugars or some such - deal with it below */
		/* by checking the backbone atoms */
		///if ( !sequence[resctr].res_type_short ) return 1;
	   
	    } else {
		atomctr ++;
		sequence[resctr].no_atoms = atomctr + 1;
		if ( atomctr >= MAX_NO_ATOMS ) {
		    fprintf ( stderr,
			      "Error parsing pdb: I thought every aa has < %d atoms.\n",
			      MAX_NO_ATOMS );
		    return ERR_MAX_ATOMS;
		}
	    }
	    /* read in atom info */
	    
	    auxptr = line+ PDB_ATOM_ATOM_NAME;
	    memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	    /* skip initial blanks*/
	    ctr  = 0;
	    while ( !(isalpha (*(auxptr + ctr))) &&
		    (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	    /* copy alphanum info */
	    nonblank = 0;
	    while (  isalnum (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
		tmp[nonblank] =  *(auxptr +ctr);
		nonblank ++;
		ctr++;
	    }

	    /* have we already seen this atom type by any chance? */
	    tmp[nonblank] = '_';
	    if ( strstr (atomtypes_read_in, tmp) ) {
		/* ahould I check for an alt location code, or just move on? */
		//printf ( " %s >> %s  //// %s\n", sequence[resctr].pdb_id, atomtypes_read_in, tmp);
		continue;
	    } else {
		sprintf (atomtypes_read_in, "%s%s", atomtypes_read_in, tmp);
	    }
	    tmp[nonblank] = '\0';
	    
	    strncpy ( sequence[resctr].atom[atomctr].type, tmp, PDB_ATOM_ATOM_NAME_LEN );

	    /* is this a backbone atom?*/
	    sequence[resctr].atom[atomctr].backbone = 0;
	    if ( nonblank == 1) {
		  sequence[resctr].atom[atomctr].backbone =
		      !(  strcmp ( tmp, "N") && strcmp ( tmp, "C") && strcmp ( tmp, "O")  );
	    } else if (  nonblank == 2) {
		if (  ! strcmp ( tmp, "CA" )) {
		    sequence[resctr].atom[atomctr].backbone = 1;
		    sequence[resctr].Ca = sequence[resctr].atom+atomctr;
		}  else  {
		    sequence[resctr].atom[atomctr].backbone = 0;
		}
	    }
	    /* check if this is Ca trace */
	    if ( strcmp ( tmp, "CA" ) ) ca_trace = 0;
	    
	    strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	    tmp[PDB_ATOM_X_LEN] = '\0';
	    sequence[resctr].atom[atomctr].x=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	    tmp[PDB_ATOM_Y_LEN] = '\0';
	    sequence[resctr].atom[atomctr].y=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	    tmp[PDB_ATOM_Z_LEN] = '\0';
	    sequence[resctr].atom[atomctr].z=atof(tmp);
	   
	}
    }

    if ( ca_trace) return ERR_SSE_NONE|ERR_CA_TRACE;
    
    /* clean PDB id tags from spaces */
    for (resctr=0; resctr < no_res; resctr ++ ) {
	retval = string_clean (sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN+1);
	if ( retval ) {
	    fprintf (stderr, "Error in read_pdb(): empty id string for residue with sequential no %d.\n", resctr);
	    return ERR_NONSENSE;
	}
    }

    /* store the sequence and its length */
    protein->sequence = sequence;
    protein->length   = no_res;

    return 0;
}


char single_letter ( char code[]){

    switch ( code[0] ) {
    case 'A':
	switch ( code [1]) {
	case 'L':
	    return 'A';
	    break;
	case 'R':
	    return 'R';
	    break;
	case 'S':
	    switch ( code[2] ) {
	    case 'N':
		return 'N';
		break;
	    case 'P':
		return  'D';
		break;
	    }
	    break;
	}
	break;
    case 'C':
	return 'C'; 
	break;
    case 'G':
	/* the second letter is always L */ 
	switch ( code[2] ) {
	case 'U':
	    return 'E';
	    break;
	case 'N':
	    return  'Q';
	    break;
	case 'Y':
	    return 'G';
	    break;
	}
	break;
    case 'H':
	return  'H';
	break;
    case 'I':
	return  'I';
	break;
    case 'L':
	switch ( code [1]) {
	case 'E':
	    return 'L';
	    break;
	case 'Y':
	    return 'K';
	    break;
	}
	break;
    case 'M':
	return 'M';
	break;
    case 'P':
	switch ( code [1]) {
	case 'H':
	    return 'F';
	    break;
	case 'R':
	    return 'P';
	    break;
	}
	break;
    case 'S':
	return 'S';
	break;
    case 'T':
	switch ( code [1]) {
	case 'H':
	    return 'T';
	    break;
	case 'R':
	    return 'W';
	    break;
	case 'Y':
	    return 'Y';
	    break;
	}
	break;
    case 'V':
	return 'V';
	break;
	
    }


    //fprintf (stdout, "Unrecognized amino acid code: %s.\n", code);
    return 0;
}

/**********************************************************/
/* get rid of spaces in a string */
int  string_clean ( char* string, int length) {
    int ctr;
    for (ctr = 0; ctr < length; ctr ++) {
	if ( isspace (string[ctr]) ) string[ctr] = '\0';
    }
    ctr=0;
    while ( !string[ctr] && ctr < length) ctr++;
    
    if ( ctr == length ) return 1; /* empty string */
    
    if ( ctr ) {
	memmove (string, string+ctr, length-ctr);
	memset ( string+length-1-ctr, 0, ctr);
    }

    return 0;
}


/***********************************************************************/
int  struct_almt_mapping (Protein * protein, Alignment * alignment, int sequence_id){
    
    int prot_pos, almt_pos;
    char * struct_seq;
    Residue * prot_seq;
    int * prot2almt = protein->prot2almt;
    int * almt2prot = protein->almt2prot;
    
    /* locate query in the alignment */
    struct_seq = alignment->sequence[sequence_id];

    /*compare */
    prot_pos = 0;
    prot_seq = protein->sequence;
    for (almt_pos=0; almt_pos < alignment->length; almt_pos++ ) {
	
	if ( struct_seq [almt_pos] == '.' ) {
	    almt2prot [almt_pos] = -1;
	} else {
	    if ( prot_seq[prot_pos].res_type_short ==  struct_seq [almt_pos] ) {
		prot2almt[prot_pos] = almt_pos;
		almt2prot[almt_pos] = prot_pos;
	    } else {
		fprintf (stderr, "Structure/alignment mismatch,\n");
		fprintf (stderr, "\t structure: pdbid %s  value %c \n",
			 prot_seq[prot_pos].pdb_id,  prot_seq[prot_pos].res_type_short);
		fprintf (stderr, "\t alignment: pos %d  value  %c \n",   almt_pos,  struct_seq[almt_pos]);
		return 1;
	    }
	    prot_pos++;
	}
    }
    return 0;
}
