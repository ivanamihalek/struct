# include "struct.h"

int results_out (Descr *tgt_descr, Protein *tgt_structure, Representation * tgt_rep,
		 Descr *qry_descr, Protein *qry_structure, Representation * qry_rep,
		 List_of_maps *list, FILE * digest) {

    
    if (options.verbose) {
	write_maps (stdout, tgt_descr, qry_descr, list);
    }

    if (options.postprocess) {
	write_tfmd_pdb  (tgt_structure, list, tgt_descr, qry_descr);
	write_alignment (tgt_structure, qry_structure, list);
    }
    
    write_digest (qry_descr, tgt_descr, qry_rep, tgt_rep, list, digest);               

    return 0;
}

