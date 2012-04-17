# include "struct.h"


int protein_shutdown (Protein * protein) {

    if ( ! protein) return 0;
    
    protein->length = 0;
    if (protein->sequence) free(protein->sequence);
    
    if (protein->no_helices) protein->no_helices = 0;
    
    if (protein->helix) free (protein->helix);
    
    if (protein->no_strands) protein->no_strands = 0;
     
    if (protein->strand) free (protein->strand);
    if (protein->sse_sequence) free (protein->sse_sequence);

    
    return 0;
}
