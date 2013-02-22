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



int process_input_instructions (int argc, char *argv[],
				int * tgt_input_type_ptr, char * tgt_chain_ptr, Descr * tgt_descr, FILE ** tgt_fptr_ptr,
				int * qry_input_type_ptr, char * qry_chain_ptr, Descr * qry_descr, FILE ** qry_fptr_ptr) {


    char tgt_chain = '\0', qry_chain = '\0';
    int  tgt_input_type = 0, qry_input_type = 0;

    FILE *qry_fptr    = NULL, *tgt_fptr = NULL;

    
    char  *cmd_filename = NULL;

    int parse_cmd_line (int argc, char * argv[],  char * tgt_chain_ptr,
			char * qry_chain_ptr, char **cmd_filename_ptr);
    int read_cmd_file (char *filename);

    
    /* process cmd line input */
    if (parse_cmd_line (argc, argv, &tgt_chain,
		        &qry_chain, &cmd_filename)) return 1;

    /* process the command (parameters) file, if provided */
    if (cmd_filename && read_cmd_file(cmd_filename))  return 1;

   
    /* check if the tgt file is  present and readable; open               */
    if ( ! (tgt_fptr = efopen(options.tgt_filename, "r"))) return 1;
    
    
    /* figure out whether we have a pdb or db input:                      */
    tgt_input_type = check_input_type (tgt_fptr);
    if ( tgt_input_type != PDB && tgt_input_type != DB ) {
	fprintf ( stderr, "Unrecognized file type: %s.\n", argv[1]);
	return 1;
    }

    /* for testing purposes we might have both: */
    if ( options.tgt_db) tgt_descr->db_file  =  options.tgt_db;
    if ( options.qry_db) qry_descr->db_file  =  options.qry_db;
    
    
    /*do something about the names for the output:                        */
    if ( tgt_input_type==PDB) {
	improvize_name (options.tgt_filename, tgt_chain, tgt_descr->name);
    }
    
    /**********************************************************************/
    /* the same for the qry file, but may not be necessary if we are      */
    /* preprocessing only                                                 */
    if ( options.qry_filename) {
	if ( ! (qry_fptr = efopen(options.qry_filename, "r"))) return 1;
	
	qry_input_type = check_input_type (qry_fptr);
	if ( qry_input_type != PDB  &&  qry_input_type != DB ) {
	    fprintf ( stderr, "Unrecognized file type: %s.\n", argv[2]);
	    exit (1);
	}
	if ( qry_input_type==PDB) {
	    improvize_name (options.qry_filename, qry_chain, qry_descr->name);
	}
	if (tgt_input_type != PDB  ||  qry_input_type != PDB)  {
	    options.postprocess        = 0;
	    options.number_maps_out = 1;
	}
	
    }

    *tgt_input_type_ptr = tgt_input_type;
    *qry_input_type_ptr = qry_input_type;

    *tgt_chain_ptr = tgt_chain;
    *qry_chain_ptr = qry_chain;

    *tgt_fptr_ptr = tgt_fptr;
    *qry_fptr_ptr  = qry_fptr;

    return 0;

}


/**************************************************************************************/
/**************************************************************************************/
int parse_cmd_line (int argc, char * argv[], char * tgt_chain_ptr,
		    char * qry_chain_ptr, char **cmd_filename_ptr) {

    int argi;

    *tgt_chain_ptr = '\0';
    *qry_chain_ptr = '\0';
   
    *cmd_filename_ptr = NULL;

    argi = 1;

    while (argi<argc) {
	
	if ( argv[argi][0] != '-' ) {
	    fprintf (stderr, "An option should be preceded by a flag: %s\n",  argv[argi]);
	    return 1;
	} else if ( ! strncmp (argv[argi], "-gpu", 4)) {
	    options.gpu = 1;
	    argi += 1;

	} else if ( ! strcmp (argv[argi], "-no_bb")) {
	    options.postprocess  = 0;
	    options.print_header = 0;
	    argi += 1;
	    
	} else if ( ! strncmp (argv[argi], "-v", 2)) {
	    options.verbose = 1;
	    argi += 1;
	    
	} else if (  argi+1 >= argc ) {
	    fprintf (stderr, "Option %s should be followed by an argument\n",  argv[argi]);
	    return 1;
	    
	} else if ( ! strcmp (argv[argi], "-in") ||  ! strncmp (argv[argi], "-in1", 4)) {
	    options.tgt_filename = argv[argi+1];
	    argi += 2;
	    
	} else if ( ! strncmp (argv[argi], "-in2", 4)) {
	    options.qry_filename = argv[argi+1];
	    argi += 2;

	} else if ( ! strncmp (argv[argi], "-c1", 3)) {
	    *tgt_chain_ptr = argv[argi+1][0];
	    argi += 2;

	} else if ( ! strncmp (argv[argi], "-c2", 3)) {
	    *qry_chain_ptr = argv[argi+1][0];
	    argi += 2;

	} else if ( ! strncmp (argv[argi], "-p", 2)) {
	    *cmd_filename_ptr = argv[argi+1];
	    argi += 2;

	/* this is upposed to be for testing purposes only: an externally provided db file */
	} else if ( ! strncmp (argv[argi], "-db1", 4)) {
	    options.tgt_db = argv[argi+1];
	    argi += 2;

	} else if ( ! strncmp (argv[argi], "-db2", 4)) {
	    options.qry_db = argv[argi+1];
	    argi += 2;	    
	    
	} else {
	    fprintf (stderr, "Unrecognized option: %s\n",  argv[argi]);
	    return 1;
	}
	

    }

    if  (!options.qry_filename)
	options.preproc_only = 1;/*automatically assume preprocessing only */
    
    
    return 0;

}



