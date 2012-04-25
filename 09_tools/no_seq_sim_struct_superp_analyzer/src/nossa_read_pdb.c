# include "nossa.h"

int read_pdb ( char * pdbname, Protein * protein, char chain) {

    Residue * sequence;
    FILE * fptr = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+2]; /* res name: 4 digits + insertion code + \0 */
    char oldrestype [PDB_ATOM_RES_NAME_LEN+2];
    char tmp[PDB_ATOM_X_LEN+1], *auxptr;
    int atomctr, resctr, no_res, ctr, nonblank;

    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return 1;
    }

    /* count residues */
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    resctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( chain && line[PDB_ATOM_CHAINID] != chain ) continue;
	if ( ! strncmp(line,"TER", 3) ||  ! strncmp(line,"END", 3) ) break;
	if( ! strncmp(line,"ATOM", 4)){
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
		/* printf ( "New residue number:  %s \n", oldresno); */
		resctr ++;
	    }
	}
    }
    no_res = resctr;
    /* printf ("no residues: %d\n", no_res); */

    /* allocate space */
    sequence = NULL;
    if ( ! (sequence = emalloc ( no_res*sizeof (Residue)))) exit (1);

    /* read in the atom */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    /*  tyring to account for the insertion code */
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    resctr= -1;
    atomctr = 0;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( chain && line[PDB_ATOM_CHAINID] != chain ) continue;
	if ( ! strncmp(line,"TER", 3) ||  ! strncmp(line,"END", 3) ) break;
	if( ! strncmp(line,"ATOM", 4) ){
	   /* if it's a hydrogen - skip */
	    if ( line[PDB_ATOM_ATOM_NAME] == 'H'
		 ||  line[PDB_ATOM_ATOM_NAME+1] == 'H') continue;
	    /* adjust the counters */ 
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		strncpy (oldrestype, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
		oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
		resctr ++;
		atomctr = 0;
		
		sequence[resctr].no_atoms = 1;
		strncpy ( sequence[resctr].pdb_id, oldresno, PDB_ATOM_RES_NO_LEN+2);
		sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN+1]   = '\0';
		
		strncpy ( sequence[resctr].res_type, oldrestype, PDB_ATOM_RES_NAME_LEN+1);
		sequence[resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';
		sequence[resctr].res_type_short  = single_letter ( sequence[resctr].res_type );
		if ( !sequence[resctr].res_type_short ) return 1;
	   
	    } else {
		atomctr ++;
		sequence[resctr].no_atoms = atomctr + 1;
		if ( atomctr >= MAX_NO_ATOMS ) {
		    fprintf ( stderr, "Error: I thought every aa has < %d atoms.\n",
			      MAX_NO_ATOMS );
		    return 1;
		}
	    }
	    /* read in atom info */
	    
	    auxptr = line+ PDB_ATOM_ATOM_NAME;
	    memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	    /* skip initial blanks*/
	    ctr  = 0;
	    while ( !(isalpha (*(auxptr + ctr))) &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	    /* copy alphanum info */
	    nonblank = 0;
	    while (  isalpha (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
		tmp[nonblank] =  *(auxptr +ctr);
		nonblank ++;
		ctr++;
	    }

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
	    /* printf ( " %4d %4d %4s is backbone: %1d \n", resctr, atomctr, */
		     /* sequence[resctr].atom[atomctr].type, sequence[resctr].atom[atomctr].backbone); */
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
 
    /* close file */
    fclose (fptr);

    /* clean PDB id tags from spaces */
    for (resctr=0; resctr < no_res; resctr ++ ) {
	string_clean (sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN+1);
    }
    
    /*return values: */
    protein->sequence= sequence;
    protein->length  = no_res;

    return 0;
}

