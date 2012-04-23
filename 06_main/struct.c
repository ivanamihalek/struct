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

Options options;

/*****************************************************************************/
/*****************************************************************************/
int main ( int argc, char * argv[]) {


    char tgt_chain = '\0', qry_chain = '\0';
    char *tgt_filename = NULL, *qry_filename = NULL, *cmd_filename = NULL;
    int retval, qry_done, tgt_done;
    int db_ctr, db_effective_ctr;
    int tgt_input_type = 0, qry_input_type = 0;
    clock_t CPU_time_begin, CPU_time_end;
    FILE *qry_fptr    = NULL, *tgt_fptr = NULL, *digest = NULL;
    Protein qry_structure = {0};
    Protein tgt_structure = {0};
    Descr qry_descr   = {{0}};
    Descr tgt_descr   = {{0}};
    Representation qry_rep = {0};
    Representation tgt_rep = {0};
    
    List_of_maps list = {NULL};
    List_of_maps list_sequential = {NULL};
    List_of_maps list_out_of_order = {NULL};
    Score score;
    Score score_sequential; 
    Score score_out_of_order;
    
    int compare_reduced_reps (Representation *rep1, Representation *rep2, List_of_maps *list, Score *score);   
    int parse_cmd_line (int argc, char * argv[], char **tgt_filename_ptr, char * tgt_chain_ptr,
			char **qry_filename_ptr, char * qry_chain_ptr, char **cmd_filename_ptr);
    int read_cmd_file (char *filename);
    int set_default_options ();
    
    if ( argc < 2 ) {
	fprintf ( stderr, "Usage: %s -in1 <pdb/db tgt file> [-c1 <tgt chain>] "
		  "[ -in2 <pdb/db qry file>] [ -c2 <qry chain>] [ -p <parameter file>].\n",
		  argv[0]);
	exit (1);
    }

    /* set defaults: */
    set_default_options ();

    /* process cmd line input */
    if (parse_cmd_line (argc, argv, &tgt_filename, &tgt_chain,
			&qry_filename, &qry_chain, &cmd_filename)) return 1;

    /* process the command (parameters) file, if provided */
    if (cmd_filename && read_cmd_file(cmd_filename))  return 1;

   
    /* check if the tgt file is  present and readable; open               */
    if ( ! (tgt_fptr = efopen(tgt_filename, "r"))) return 1;
    /* figure out whether we have a pdb or db input:                      */
    tgt_input_type = check_input_type (tgt_fptr);
    if ( tgt_input_type != PDB && tgt_input_type != DB ) {
	fprintf ( stderr, "Unrecognized file type: %s.\n", argv[1]);
	exit (1);
    }
    /*do something about the names for the output:                        */
    if ( tgt_input_type==PDB) {
	improvize_name (tgt_filename, tgt_chain, tgt_descr.name);
    }
    
    /**********************************************************************/
    /* the same for the qry file, but may not be necessary if we are      */
    /* preprocessing only                                                 */
    if ( qry_filename) {
	if ( ! (qry_fptr = efopen(qry_filename, "r"))) return 1;
	qry_input_type = check_input_type (qry_fptr);
	if ( qry_input_type != PDB  &&  qry_input_type != DB ) {
	    fprintf ( stderr, "Unrecognized file type: %s.\n", argv[2]);
	    exit (1);
	}
	if ( qry_input_type==PDB) {
	    improvize_name (qry_filename, qry_chain, qry_descr.name);
	}
	if (options.postprocess) {
	    if (tgt_input_type != PDB  ||  qry_input_type != PDB) {
		fprintf ( stderr, "Both input files must be PDB to do the postprocessing.\n");
		exit (1);
	    }
	}
    }

    if (options.preproc_only) {
	/***********************************************************************/
	/***********************************************************************/
	/***********************************************************************/
	/* preprocessing only :                                                */
	tgt_done = 0;
	db_ctr   = 0;
	while ( ! tgt_done) {
	    db_ctr++;
	    retval = get_next_descr (tgt_input_type, tgt_fptr, tgt_chain, &tgt_structure, &tgt_descr);
	    if ( retval != 0  ) { /* I don't know how to recover if this is a concat of PDB files      */
				  /* actually if I want to concatenate I have another problem: rewind  */
		                  /* so PDB input should be completley                                 */
				  /* rewritten to acommodate this possibility                          */
		tgt_done = 1;
		continue;
	    }
	    descr_out (NULL, &tgt_descr);
 	}
   
    } else {	
	/**********************************************************************/
	/* read in the table of integral values                               */
	/* the array int_table in struct_table.c                              */
	if ( read_integral_table (options.path) ) {
	    fprintf (stderr, "In data file  %s.\n\n", options.path);
	    exit (1);
	}
	set_up_exp_table ();
   
	/***********************************************************************/
	/* initialize the digest file: concise report on the comparison scores */
	/*   for each pair of the structures we are looking at                 */
	init_digest (&qry_descr, &tgt_descr, &digest);
 	/***********************************************************************/
	/***********************************************************************/
	/***********************************************************************/
	/* compare pairs from tgt and qry lists :                              */
	list_alloc (&list, INIT_ALLOC_N, INIT_ALLOC_N);
	qry_done = 0;
	retval = -1;
	db_effective_ctr = 0;
	CPU_time_begin = clock();   
	while ( ! qry_done) {
            
	    retval = get_next_descr (qry_input_type, qry_fptr, qry_chain, &qry_structure, &qry_descr);
	    if ( retval == 1 ) {
		continue;
	    } else if ( retval == -1 ) {
		qry_done = 1;
		continue;
	    }

	    /*******************************/
	    /* loop over target  database :*/
	    rewind (tgt_fptr);
	    tgt_done = 0;
	    db_ctr   = 0;
	    db_effective_ctr = 0;
	    CPU_time_begin = clock();
	    retval = -1;
    
	    while ( ! tgt_done) {
		db_ctr++;
		retval = get_next_descr (tgt_input_type, tgt_fptr, tgt_chain, &tgt_structure, &tgt_descr);
		if ( retval == 1 ) {
		    continue;
		} else if ( retval == -1 ) {
		    tgt_done = 1;
		    continue;
		} 

		int match_found = 0;
		/* min number of elements */
		int helix_overlap =
		    (qry_descr.no_of_helices < tgt_descr.no_of_helices) ?
		    qry_descr.no_of_helices : tgt_descr.no_of_helices;
		int strand_overlap =
		    (qry_descr.no_of_strands < tgt_descr.no_of_strands) ?
		    qry_descr.no_of_strands : tgt_descr.no_of_strands;

		memset (&score, 0, sizeof(Score));
		
		if ( helix_overlap + strand_overlap >= options.min_no_SSEs) {

		    rep_initialize (&tgt_rep, &tgt_descr);
		    rep_initialize (&qry_rep, &qry_descr);
    
		    /*************************************************************/
		    /*************************************************************/
		    /*  here is the core: comparison of reduced representations  */
                    int retval1 = 0, retval2 = 0;
                    
                    switch (options.search_algorithm){
                        case sequential:
                            options.current_algorithm = sequential;    
                            retval1 = compare_reduced_reps ( &tgt_rep, &qry_rep, &list_sequential, &score_sequential);
                            list = list_sequential;
                            score = score_sequential;
                            break;
                        case out_of_order:
                            options.current_algorithm = out_of_order;    
                            retval1 = compare_reduced_reps ( &tgt_rep, &qry_rep, &list_out_of_order, &score_out_of_order);
                            list = list_out_of_order;
                            score = score_out_of_order;
                            break;
                        case both:
                            options.current_algorithm = sequential;    
                            retval1 = compare_reduced_reps ( &tgt_rep, &qry_rep, &list, &score);
                            options.current_algorithm = out_of_order;    
                            retval1 = compare_reduced_reps ( &tgt_rep, &qry_rep, &list_out_of_order, &score_out_of_order);
                            if (score_sequential.total_assigned_score > score_out_of_order.total_assigned_score) {
                                list = list_sequential;
                                score = score_sequential;
                            } else {
                                list = list_out_of_order;
                                score = score_out_of_order;
                            }
                    }
                    
		    if (retval1 || retval2) { /* this might be printf (rather than fprintf)
				     bcs perl has a problem intercepting stderr */
			printf (" error comparing   db:%s  query:%s \n",
				tgt_descr.name, qry_descr.name);
			exit (1);
		    }
		    db_effective_ctr ++;
                    write_digest(&qry_descr, &tgt_descr, digest, &score);               
		   
		   
		    if (options.verbose) {
			write_maps (stdout, &tgt_descr, &qry_descr, &list);
		    }
		    match_found = (list.map_best[0] >= 0);

		    if  (match_found) { /* otherwise we had no match */
			if (options.postprocess) {
			    retval = align_backbone  (&tgt_descr, &tgt_structure, &tgt_rep,
					 &qry_descr, &qry_structure, &qry_rep, &list, &score);
			    if (  retval) {
				printf (" error doing bb alignment   db:%s  query:%s \n",
					tgt_descr.name, qry_descr.name);
				exit (retval);
			    }
			    write_tfmd_pdb  (&tgt_structure, &list, &tgt_descr, &qry_descr);
			    write_alignment (&tgt_structure, &qry_structure, &list);
			}
			printf ("match found for db:%s  query:%s \n",
				tgt_descr.name, qry_descr.name);
		    } else {
			printf ("no match for db:%s  query:%s \n",
				tgt_descr.name, qry_descr.name);
		    }
		    
		    rep_shutdown (&tgt_rep);
		    rep_shutdown (&qry_rep);
    
		} else if (options.report_no_sse_overlap) {
		    /* write all zeros to the digest file  */
		    write_digest(&qry_descr, &tgt_descr, digest, &score);
		    printf ("no common SSEs for db:%s  query:%s \n",
				tgt_descr.name, qry_descr.name);
		}
	    }
	}

	CPU_time_end = clock();
	close_digest(CPU_time_begin, CPU_time_end, digest);
 
	if (options.verbose ) {
	    printf ("\n\nlooked at %d db entries.\n", db_effective_ctr);
	    printf ("CPU:  %10.3lf s\n", (double)(CPU_time_end-CPU_time_begin)/CLOCKS_PER_SEC );
	    printf ("the output written to %s.\n\n", options.outname);
	}
	
	list_shutdown (&list); /* defined in struct_map */
    }

    descr_shutdown (&qry_descr);
    descr_shutdown (&tgt_descr);
    
    protein_shutdown (&qry_structure);
    protein_shutdown (&tgt_structure);

    if (qry_fptr) fclose (qry_fptr);  
    if (tgt_fptr) fclose (tgt_fptr);
    
    return 0;
    
}

/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
int parse_cmd_line (int argc, char * argv[], char **tgt_filename_ptr, char * tgt_chain_ptr,
		    char **qry_filename_ptr, char * qry_chain_ptr, char **cmd_filename_ptr) {

    int argi;

    *tgt_filename_ptr =  NULL;
    *qry_filename_ptr =  NULL;

    *tgt_chain_ptr = '\0';
    *qry_chain_ptr = '\0';
   
    *cmd_filename_ptr = NULL;


    for (argi=1; argi<argc; argi+=2) {
	
	if ( argv[argi][0] != '-' ) {
	    fprintf (stderr, "An option should be preceded by a flag: %s\n",  argv[argi]);
	    return 1;
	} else if (  argi+1 >= argc ) {
	    fprintf (stderr, "Option %s should be followed by an argument\n",  argv[argi]);
	    return 1;
	} else if ( ! strncmp (argv[argi], "-in1", 4)) {
	    *tgt_filename_ptr = argv[argi+1];
	} else if ( ! strncmp (argv[argi], "-in2", 4)) {
	    *qry_filename_ptr = argv[argi+1];
	} else if ( ! strncmp (argv[argi], "-c1", 3)) {
	    *tgt_chain_ptr = argv[argi+1][0];
	} else if ( ! strncmp (argv[argi], "-c2", 3)) {
	    *qry_chain_ptr = argv[argi+1][0];
	} else if ( ! strncmp (argv[argi], "-p", 2)) {
	    *cmd_filename_ptr = argv[argi+1];
	} else {
	    fprintf (stderr, "Unrecognized option: %s\n",  argv[argi]);
	    return 1;
	}
	

    }

    if  (!*qry_filename_ptr)
	options.preproc_only = 1;/*automatically assume preprocessing only */
    
    
    return 0;

}



/**************************************************************************/

int set_default_options () {
    /* set the default options */
    memset (&options, 0, sizeof(Options) );


    options.min_no_SSEs = 4;

    options.merge_cosine = 1.1; /*min cos angle between two SSEs to be represented
				  by (merged into) the same vector */
    
    options.alpha = 0.0;  /* gaussian width for the scoring fn F */
                          /* note that it is set from the input table */
                          /* rather than from the cmd file        */
    options.tol           /* tol in F for the claim that the query and  */
	                  /* the target are one and the same structure  */
	= 0.1;            
    options.z_max_store   /* z-score cutoff to store a map */
	= -2.0;
    options.z_max_out     /* z-score cutoff for map output  */
	= -2.0;
	
    options.F_guess_max   /* max  value of F for an intial guess    */
	= -1.5;
    options.F_eff_max     /* max  value of F after gradient descent */
	= -3.0;
    options.z_max_corr    /* max value of z for which two maps are */
			  /* considered correlated */
	= -2;
    options.z_min_compl   /* min value of z for which two maps are */
			  /* considered complementary */
	=  2;
    options.grad_step_size/* step size in gradient descent */
	= 0.0001;
    options.grad_stop_tol /* stopping precision in gradient descent */
	= 1.e-4;
    
    options.far_far_away  /* nonsense distance for the pairwise alignment*/
	= -10;
    options.gap_open      /* gap penalties for the pairwise alignment */
	= 0.0;
    options.gap_extend
	= 0.0;
    options.endgap
	= 0.0;
    options.threshold_distance
	= 15.0;
    options.far_away_cosine /* minimum cosine for F_effective estimate*/
        = 0.8;
    options.grid_size     /* minimal number of points for the sphere grid */
	= 400;
    options.number_maps_cpl /* number of top scoring maps among which to   */
                            /* look for a complement */ 
	= 0;
    options.number_maps_out /* number of top scoring maps to output     */
	= 0;
    options.grad_max_step   /* max number of steps in  gradient descent */
	= 100;
    options.exp_table_size  /* size of the lookup table for the exp funcion */
	= TABLE_SIZE;       /* note: changing exp table size not implemented */

    options.exhaustive     = 0; /* try all triples instead of consecutive only */
    options.smith_waterman = 1;

    options.search_algorithm = sequential;
    
    
    options.verbose        = 0;
    options.print_header   = 0;
    options.report_no_sse_overlap = 0;
    
    /* penalizing the length mismatch */
    options.use_length = 0;
    options.H_length_mismatch_tol = 10.0;
    options.S_length_mismatch_tol = 5.0;
    
    /* path to the integral table */
    memset (options.path, 0, BUFFLEN);
    sprintf (options.path,  "%s",  INTEGRAL_TABLE);

    return 0;
}


