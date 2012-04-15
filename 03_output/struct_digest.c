# include "struct.h"

int init_digest (Descr *qry_descr, Descr *tgt_descr, FILE ** digest_ptr) {


    FILE *digest = *digest_ptr;
    
    if (!digest) { /*open new one */
	char outname[MEDSTRING] = {'\0'};
	if (!options.outname[0] ) {
	    if ( qry_descr->name[0] &&  tgt_descr->name[0] ) {
		sprintf (outname, "%s_%s.struct_out", qry_descr->name, tgt_descr->name);
	    } else {
		sprintf (outname, "digest.struct_out");
	    }
	} else {
	    sprintf (outname, "%s.struct_out", options.outname);
	}
	digest  = efopen (outname, "w");
    }
    
    if ( !digest) return 1;

    *digest_ptr = digest;
    
    if ( options.print_header) {
	fprintf ( digest,"%% columns: \n");
	fprintf ( digest,"%% query, target: structure names\n");
	fprintf ( digest,"%% geom_z:  z score for the orientational match \n");
	fprintf ( digest,"%% <dL>:    average length mismatch for matched SSEs \n");
	fprintf ( digest,"%% T:       total score assigned to matched SSEs \n");
	fprintf ( digest,"%% frac:    T divided by the number of matched SSEs \n");
	fprintf ( digest,"%% GC_rmsd: RMSD btw geometric centers "
		  "of matched SSEs (before postprocessing) \n");
	fprintf ( digest,"%% A:       (after postprocessing) the alignment score \n");
	fprintf ( digest,"%% aln_L:   (after postprocessing) the alignment length \n\n");
	fprintf ( digest,"%% %6s%6s %6s %6s  %6s %6s %6s %6s %6s %6s \n",
		  "query ", "target ", "geom_z", "<dL>", "  T  ", "frac",
		  "GC_rmsd", "rmsd  ", "A  ", "aln_L  " );	
    }

    return 0;

}


int write_digest(Descr *qry_descr, Descr *tgt_descr, FILE * digest, Score * score) {
    
    fprintf ( digest,
	      "%6s %6s %8.3lf %6.3lf %6.2lf %6.2lf %6.3lf %6.3lf %6.3lf %4d ",
	      qry_descr->name,
	      tgt_descr->name,
			      
	      score->z_score,
	      score->avg_length_mismatch,
			      
	      score->total_assigned_score,
	      score->fraction_assigned,
			  
	      score->rmsd,
			      
	      score->res_rmsd,
	      score->res_almt_score,
	      score->res_almt_length);
    fprintf (digest, "\n");
    fflush  (digest);

    return 0;

}


int close_digest (clock_t CPU_time_begin, clock_t CPU_time_end, FILE *digest){

    fprintf (digest, "done   CPU:  %10.3lf s\n", (double)(CPU_time_end-CPU_time_begin)/CLOCKS_PER_SEC );
    fflush  (digest);
    fclose(digest);

    return 0;
}
