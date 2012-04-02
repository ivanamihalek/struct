# include "struct.h"

# define MISSING_DATA  -1

int structure2seq (Protein *protein);

int main ( int argc, char * argv[]) {

    int retval;
    int seq_length = 0;
    Protein protein;
   
    if ( argc < 3 ) {
	fprintf ( stderr,
		  "Usage: %s <pdb file> [<chain>].\n",
		  argv[0]);
	exit (1);
    }

    retval =  read_pdb(argv[1], argv[2][0], &protein);
    if ( retval ) {
	fprintf ( stderr, "Error reading %s, retval %d.\n",
		  argv[1], retval);
	exit (1);
    }


    if ( structure2seq (&protein)) {
	fprintf ( stderr, "Error  creating motif sequence.\n");
	exit (1);
    }
    

    int  resctr;

    for (resctr=0; resctr<protein.length; resctr++) {

	char sse_code = 'C';
	if ( protein.sequence[resctr].belongs_to_helix) {
	    sse_code = 'H';
	} else if ( protein.sequence[resctr].belongs_to_strand) {
	    sse_code = 'S';
	}
	
	printf ("  %s  %c\n", protein.sequence[resctr].pdb_id, sse_code);

	
    }
    
    
    return 0;
}




/*************************************************************/
/* in this first version, the sequence is a sequence of distances
   to the first three upstream neighbors */
int structure2seq (Protein *protein) {

    int resctr1, resctr2, ctr;
    double x1, y1, z1;
    double x2, y2, z2;
    double aux, d;
    Residue *res1, *res2;

 
    ctr = 0;
    for (resctr1=0; resctr1<protein->length; resctr1++) {
	
	res1 = protein->sequence + resctr1;

	if ( ! res1->Ca) {
	    /* this (no C-alpha) may happen when the protein is modified,
	       either posttranslationally or by the crystallographer */
	    for (resctr2=resctr1+2; resctr2<=resctr1+4 ; resctr2++) {
		//motif_seq[ctr] = MISSING_DATA;
		ctr++;
	    }
	    continue;
	}
	
	x1 = res1->Ca->x;
	y1 = res1->Ca->y;
	z1 = res1->Ca->z;
	
	for (resctr2=resctr1+1; resctr2<protein->length ; resctr2++) {
	    
	    res2 = protein->sequence + resctr2;

	    if ( ! res2->Ca) {
		//motif_seq[ctr] = MISSING_DATA;
		
	    } else {
		x2 = res2->Ca->x;
		y2 = res2->Ca->y;
		z2 = res2->Ca->z;

		d = 0;
		aux = x1-x2; d += aux*aux;
		aux = y1-y2; d += aux*aux;
		aux = z1-z2; d += aux*aux;
	    
	    }
	    ctr++;
	}
	
    }
    
    /* fill protein->sequence[resctr].belongs_to_helix or
       protein->sequence[resctr].belongs_to_strand */

    
    
    return 0;
}
