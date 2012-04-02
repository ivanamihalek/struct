# include "struct.h"




# define DELETE  -789


int read_pdb ( char * pdbname,  char chain, Protein * protein) {

    /* TODO for the moment we rely on PDB annotation
       to extract structural elements - that should be changed */
    Residue * sequence = NULL;
    FILE    * fptr = NULL;
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
    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return ERR_NO_FILE_OR_CHAIN;
    }


    
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
	
	if ( resctr ) {
	    if ( ! strncmp(line,"TER", 3) ||  ! strncmp(line,"END", 3) ||  line[PDB_ATOM_CHAINID] != old_chain)
		break;
	}
	if ( chain  && line[PDB_ATOM_CHAINID] != chain ) continue;
	chain_found  = 1;
	
	if( ! strncmp(line,"ATOM", 4) ||  ! strncmp(line,"HETATM", 6)){
	    
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
	fprintf (stderr, "Chain %c not found in %s\n", chain, pdbname);
	return ERR_NO_FILE_OR_CHAIN;
    }
    
    no_res = resctr;
   
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
	    if  (! strncmp(line,"TER", 3) ||  ! strncmp(line,"END", 3)  ||  line[PDB_ATOM_CHAINID] != old_chain)
	    break;
	}
	if ( chain  && line[PDB_ATOM_CHAINID] != chain ) continue;
	
	if( ! strncmp(line,"ATOM", 4) ||  ! strncmp(line,"HETATM", 6)){
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
		    fprintf (stderr, "Error reading %s: rectr:%d   no res: %d\n", pdbname, resctr, no_res);
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
		    fprintf ( stderr, "Error parsing %s: I thought every aa has < %d atoms.\n",
			      pdbname, MAX_NO_ATOMS );
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
	    fprintf (stderr, "Error in read_pdb(): empty id string for residue %d.\n", resctr);
	    return ERR_NONSENSE;
	}
    }

    /* store the sequence and its length */
    protein->sequence = sequence;
    protein->length   = no_res;

    fclose(fptr);
    return 0;
}


