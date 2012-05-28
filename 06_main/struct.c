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

Options options;

/*****************************************************************************/
/*****************************************************************************/
int main ( int argc, char * argv[]) {


    char tgt_chain = '\0', qry_chain = '\0';
    int retval, qry_done, tgt_done;
    int db_ctr, db_effective_ctr;
    int tgt_input_type = 0, qry_input_type = 0;
    clock_t CPU_time_begin, CPU_time_end, CPU_comparison_start;
    FILE *qry_fptr    = NULL, *tgt_fptr = NULL, *digest = NULL;
    Protein qry_structure = {0};
    Protein tgt_structure = {0};
    Descr qry_descr   = {{0}};
    Descr tgt_descr   = {{0}};
    Representation qry_rep = {0};
    Representation tgt_rep = {0};
    
    List_of_maps list_sequential = {NULL};
    List_of_maps list_out_of_order = {NULL};
    List_of_maps list_uniq = {NULL};
    List_of_maps *list1=NULL, *list2=NULL;
     
    int map_reduced_reps (Representation *rep1, Representation *rep2, List_of_maps *list);   
    int process_input_instructions (int argc, char *argv[],
				    int * tgt_input_type_ptr, char * tgt_chain_ptr, Descr * tgt_descr, FILE ** tgt_fptr_ptr,
				    int * qry_input_type_ptr, char * qry_chain_ptr, Descr * qry_descr, FILE ** qry_fptr_ptr);
    int set_default_options ();
    
    if ( argc < 2 ) {
	fprintf ( stderr, "Usage: %s -in1 <pdb/db tgt file> [-c1 <tgt chain>] "
		  "[ -in2 <pdb/db qry file>] [ -c2 <qry chain>] [-no_bb] [ -p <parameter file>].\n",
		  argv[0]);
	exit (1);
    }

    /* set defaults: */
    set_default_options ();

    retval = process_input_instructions(argc, argv,
					&tgt_input_type, &tgt_chain, &tgt_descr, &tgt_fptr,
					&qry_input_type, &qry_chain, &qry_descr, &qry_fptr);
    if (retval) return retval;
    
    
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
	if ( options.path[0] ) {
	    if ( read_integral_table (options.path) ) {
		fprintf (stderr, "In data file  %s.\n\n", options.path);
		exit (1);
	    }
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
	int fake;
	list_alloc (&list_sequential,   INIT_ALLOC_N, INIT_ALLOC_N, (fake=0));
        list_alloc (&list_out_of_order, INIT_ALLOC_N, INIT_ALLOC_N, (fake=0));
        list_alloc (&list_uniq,         INIT_ALLOC_N, INIT_ALLOC_N, (fake=1));
        
	qry_done = 0;
	retval = -1;
	db_effective_ctr = 0;
	CPU_time_begin = clock();
	while ( ! qry_done) {
            
	    retval = get_next_descr (qry_input_type, qry_fptr, qry_chain, &qry_structure, &qry_descr);
	    if ( retval == 1 ) {
		fprintf (stderr, "Error reading %s.\n", options.qry_filename);
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
		    fprintf (stderr, "Error reading %s.\n", options.tgt_filename);
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

		
		if ( helix_overlap + strand_overlap >= options.min_no_SSEs) {

		    
		    if (options.verbose) printf ("\n\n---------------\ncomparing   db:%s  query:%s \n",
						 tgt_descr.name, qry_descr.name);
		    CPU_comparison_start = clock();
		    rep_initialize (&tgt_rep, &tgt_descr);
		    rep_initialize (&qry_rep, &qry_descr);
   
		    /*************************************************************/
		    /*************************************************************/
		    /*  here is the core: comparison of reduced representations  */
                    int retval1 = 0, retval2 = 0;
                    
                    switch (options.search_algorithm){
			
		    case SEQUENTIAL:
			options.current_algorithm = SEQUENTIAL;    
			retval1 = map_reduced_reps (&tgt_rep, &qry_rep, &list_sequential);
			match_found = (list_sequential.no_maps_used > 0);
			list1   = &list_sequential;  list2 = NULL;
			break;
			    
		    case OUT_OF_ORDER:
			options.current_algorithm = OUT_OF_ORDER;    
			retval2 = map_reduced_reps (&tgt_rep, &qry_rep, &list_out_of_order);
			match_found =  (list_out_of_order.no_maps_used > 0);
			list1   = NULL;  list2 = &list_out_of_order;
			break;
			    
		    case BOTH:
			options.current_algorithm = SEQUENTIAL;    
			retval1 = map_reduced_reps (&tgt_rep, &qry_rep, &list_sequential);
			options.current_algorithm = OUT_OF_ORDER;    
			retval2 = map_reduced_reps (&tgt_rep, &qry_rep, &list_out_of_order);
			match_found =  (list_sequential.no_maps_used > 0 ||
					list_out_of_order.no_maps_used > 0);
			list1   = &list_sequential;  list2 = &list_out_of_order;
                    }
                    
		    if (retval1 || retval2) { /* this might be printf (rather than fprintf)
						 bcs perl has a problem intercepting stderr */
			printf (" error comparing   db:%s  query:%s \n",
				tgt_descr.name, qry_descr.name);
			exit (1);
		    }
		    db_effective_ctr ++;
		    
		    printf (" db:%s  query:%s   CPU:  %10.3lf s\n", tgt_descr.name, qry_descr.name,
		    	    (double)(clock()-CPU_comparison_start)/CLOCKS_PER_SEC );

		    if  (match_found) {

			find_uniq_maps (list1, list2, &list_uniq);
			
			if (options.postprocess) {
			    align_backbone (&tgt_descr, &tgt_structure, &tgt_rep,
					    &qry_descr, &qry_structure, &qry_rep, &list_uniq);
			    if (options.optimize)  {
				optimize_backbone_alignment (&tgt_descr, &tgt_structure, &tgt_rep,
					    &qry_descr, &qry_structure, &qry_rep, &list_uniq);
			    }
			}
			
			
			results_out (&tgt_descr, &tgt_structure, &tgt_rep,
				     &qry_descr, &qry_structure, &qry_rep,
				     &list_uniq, digest);
			
			if (options.verbose)
			    printf ("match found for db:%s  query:%s \n",
						     tgt_descr.name, qry_descr.name);
			
		    } else {
			/* write all zeros to the digest file  */
			if (options.report_no_match)
			    write_digest(&qry_descr, &tgt_descr, &qry_rep, &tgt_rep, NULL, digest);
			if (options.verbose) printf ("no match for db:%s  query:%s \n",
						     tgt_descr.name, qry_descr.name);
		    }
		    
		    rep_shutdown (&tgt_rep);
		    rep_shutdown (&qry_rep);
    
		} else if (options.report_no_sse_overlap) {
		    /* write all zeros to the digest file  */
		    write_digest(&qry_descr, &tgt_descr, &qry_rep, &tgt_rep, NULL, digest);
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
	
	list_shutdown (&list_sequential,   (fake=0));   /* defined in struct_map */
	list_shutdown (&list_out_of_order, (fake=0)); /* defined in struct_map */
	list_shutdown (&list_uniq, (fake=1)); /* defined in struct_map */
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
/**************************************************************************/

int set_default_options () {
    /* set the default options */
    memset (&options, 0, sizeof(Options) );


    options.min_no_SSEs = 3;

    options.merge_cosine = 1.1; /*min cos angle between two SSEs to be represented
				  by (merged into) the same vector */
    
    options.alpha = 0.3;  /* gaussian width for the scoring fn F */
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
    
    options.far_far_away  /* nonsense distance for the pairwise alignment something like -10*/
	= FAR_FAR_AWAY;
    options.gap_open      /* gap penalties for the pairwise alignment */
	= 0.0;
    options.gap_extend
	= 0.0;
    options.endgap
	= 0.0;
    options.threshold_distance /* for the "seed" triples in direction search */
	= 30.0;
    options.far_away_cosine /* minimum cosine for F_effective estimate*/
        = 0.8;
    options.grid_size     /* minimal number of points for the sphere grid */
	= 400;
    options.number_maps_cpl /* number of top scoring maps among which to   */
                            /* look for a complement */ 
	= 0;
    options.number_maps_out /* number of top scoring maps to output     */
	= 10;
    options.grad_max_step   /* max number of steps in  gradient descent */
	= 100;
    options.exp_table_size  /* size of the lookup table for the exp funcion */
	= TABLE_SIZE;       /* note: changing exp table size not implemented */

    options.exhaustive     = 1; /* try all triples instead of consecutive only */
    options.smith_waterman = 1;

    options.search_algorithm = SEQUENTIAL;
    

    options.postprocess    = 1;
    options.optimize       = 1;
    
    options.verbose        = 0;
    options.print_header   = 1;
    options.report_no_sse_overlap = 0; /* refers to the digest file */
    options.report_no_match = 0;      /* refers to the digest file */
    
    /* penalizing the length mismatch */
    options.use_length              =  1;
    options.H_length_mismatch_tol   = 10.0;
    options.S_length_mismatch_tol   =  5.0;
    options.avg_length_mismatch_tol =  5.0;
    
    /* path to the integral table */
    memset (options.path, 0, BUFFLEN);

    return 0;
}


