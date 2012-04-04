# include "struct.h"
# include "math.h"

# define MISSING_DATA  -1

typedef struct {    
    double dist2; // distance to Ca atom with res_id + 2
    double dist3; // distance to Ca atom with res_id + 3
    double dist4; // distance to Ca atom with res_id + 4
} Neighbors;

typedef struct  {
    char struct_type;
    Neighbors ideal_rec;
    int count_min;
    double weight_thr;
} Ideal_struct;


int determine_sec_structure(Neighbors *neighbors, Protein *protein);
int clean_short_structures(Protein * protein, Ideal_struct * ideal_helix, Ideal_struct * ideal_strand);
int calculate_neighbors_distances(Protein *protein, Neighbors * neighbors);
int structure2seq (Protein *protein);
void return_missing_data_distances(Neighbors * neighbors);
double calc_dist_res(Residue *res1, Residue * res2);
char struct_type(Residue * res);
int trace_back(int count, int count_min, int id, Protein * protein);
int is_regular_struct(Neighbors * neighbors, Ideal_struct * ideal_struct);


int main ( int argc, char * argv[]) {

    int retval;
    //int seq_length = 0;
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

void return_missing_data_distances(Neighbors * neighbors) {
    neighbors->dist2 = MISSING_DATA;
    neighbors->dist3 = MISSING_DATA;
    neighbors->dist4 = MISSING_DATA;
}

double calc_dist_res(Residue *res1, Residue *res2){
    double x1, y1, z1, x2, y2, z2;
    double distance = 0, aux;
    
    x1 = res1->Ca->x;
    y1 = res1->Ca->y;
    z1 = res1->Ca->z;
	
    x2 = res2->Ca->x;
    y2 = res2->Ca->y;
    z2 = res2->Ca->z;
		
 
    aux = x1-x2; distance += aux*aux;
    aux = y1-y2; distance += aux*aux;
    aux = z1-z2; distance += aux*aux;
    
    return sqrt(distance);   
}

int is_regular_struct(Neighbors * neighbors, Ideal_struct * ideal_struct){
    double weight = 0;
    double diff2, diff3, diff4;
    diff2 = neighbors->dist2 - ideal_struct->ideal_rec->dist2;
    diff3 = neighbors->dist3 - ideal_struct->ideal_rec->dist3;
    diff4 = neighbors->dist4 - ideal_struct->ideal_rec->dist4;
    
    weight = sqrt(diff2*diff2+diff3*diff3+diff4*diff4);
    if (weight < ideal_struct->weight_thr) return 1; 
    return 0;
}


int calculate_neighbors_distances(Protein *protein, Neighbors * neighbors) {
    
    Residue res1, res2;
    int i, j;
    int start_shift = 2;
    int end_shift = 4;
    for (i=0; i<protein->length-4; i++) {
	
        int missing_Ca = 0;
	res1 = protein->sequence[i];

	if ( ! res1.Ca) {
	    /* this (no C-alpha) may happen when the protein is modified,
	       either posttranslationally or by the crystallographer */
            return_missing_data_distances(&neighbors[i]);
	    continue;
	}
	
        for (j = i + start_shift; j <= i + end_shift; j++){
            res2 = protein->sequence[j];
            if (!res2.Ca) {
                missing_Ca = 1;
                return_missing_data_distances(&neighbors[i]);
                break;
            }
        }
        
        if (missing_Ca) continue;
        
        neighbors[i].dist2 = calc_dist_res(&res1, &protein->sequence[i+2]);
        neighbors[i].dist3 = calc_dist_res(&res1, &protein->sequence[i+3]);
        neighbors[i].dist4 = calc_dist_res(&res1, &protein->sequence[i+4]);
    }
    return 1;
}

int trace_back(int count, int count_min, int id, Protein * protein) {
    int length = protein->length;
    int i;
    
    if (count < count_min) {
        for (i = 0; i < length; ++i) {
            protein->sequence[id - i].belongs_to_helix  = 0;
            protein->sequence[id - i].belongs_to_strand = 0;
        }
    }
    return 1;
}


int clean_short_structures(Protein * protein, Ideal_struct * ideal_helix, Ideal_struct * ideal_strand) {
    char state_old = 'S', state;
    int count = 0;
    int res_id_old = 0, res_id;
    int i;
    int length = protein->length;
    
    for (i = 0; i < length; ++i) {
        state = struct_type(&protein->sequence[i]);
        switch (state) {
            case 'C':
                if (state_old == 'H') {
                    trace_back(count, ideal_helix->count_min, i -1, protein );
                } 
                else if (state_old == 'E') {
                    trace_back(count, ideal_strand->count_min, i -1, protein );
                }
                // in the case of 'C' or 'S' it is not necessary to do traceback
                count = 0;
                state_old = 'C';
                break;
                
            case 'H':
                if (state_old == 'S' || state_old == 'C') {
                    count = 1;
                }
                else if (state_old == 'E') {
                    trace_back(count, ideal_strand->count_min, i -1, protein );
                    count = 1;
                } else if (state_old == 'H') {
                    if (res_id_old + 1 == res_id){
                        count++;
                    }
                    else {
                         trace_back(count, ideal_helix->count_min, i -1, protein );
                         count = 0;
                    }
                }
                res_id_old = res_id;
                state_old = 'H';
                break;
                
            case 'E':
                if (state_old == 'S' || state_old == 'C') {
                    count = 1;
                }
                else if (state_old == 'H') {
                    trace_back(count, ideal_helix->count_min, i -1, protein );
                    count = 1;
                } else if (state_old == 'E') {
                    if (res_id_old + 1 == res_id){
                        count++;
                    }
                    else {
                         trace_back(count, ideal_strand->count_min, i -1, protein );
                         count = 0;
                    }
                }
                res_id_old = res_id;
                state_old = 'E';
                
        }
    }
    return 1;
    
}

char struct_type(Residue * res) {
    char type = 'C';
    if (res->belongs_to_helix) type = 'H';
    else if (res->alt_belongs_to_helix) type = 'E';
    
    return type;
}


int determine_sec_structure(Neighbors *neighbors, Protein *protein) {
    int length = protein->length;
    int is_struct;
    int i;
    Ideal_struct ideal_helix = {
        .struct_type = 'H',
        .ideal_rec = {5.43, 5.05, 6.20},
        .count_min = 4,
        .weight_thr = 1.5
    };
    Ideal_struct ideal_strand = {
        .struct_type = 'E',
        .ideal_rec = {6.2, 9.5, 12.4},
        .count_min = 3,
        .weight_thr = 2
    };
    
    
    for (i=0; i < length; ++i) {
        if (neighbors[i].dist2 == MISSING_DATA) continue;
        is_struct = is_regular_struct(neighbors +i, &ideal_helix);
        if (is_struct) {
            protein->sequence[i].belongs_to_helix = 1;
            continue;
        }
        is_struct = is_regular_struct(neighbors +i, &ideal_strand);
        if (is_struct) {
            protein->sequence[i].belongs_to_strand = 1;
        }
        
    }
    
    clean_short_structures(protein, &ideal_helix, &ideal_strand);
    
    return 1;
}




/*************************************************************/
/* in this first version, the sequence is a sequence of distances
   to the first three upstream neighbors */
int structure2seq (Protein *protein) {

    Neighbors* neighbors = malloc(sizeof(Neighbors)*protein->length);
    calculate_neighbors_distances(protein, neighbors);
    
    determine_sec_structure(neighbors, protein);
    free(neighbors);
    
    
    
    
    /* fill protein->sequence[resctr].belongs_to_helix or
       protein->sequence[resctr].belongs_to_strand */

    
    
    return 1;
}
