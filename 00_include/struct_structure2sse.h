# include "struct.h"

# define MISSING_DATA  -1

typedef struct {    
    double dist2; // distance to Ca atom with res_id + 2
    double dist3; // distance to Ca atom with res_id + 3
    double dist4; // distance to Ca atom with res_id + 4
} Neighbors;

typedef struct  {
    char struct_type;
    Neighbors ideal_rec;  // distances to upstream neighbors for an ideal structure
    int count_min; // min number of consecutive secondary structures
    double weight_thr; //  threshold value for the optimization function
} Ideal_struct;


int determine_sec_structure(Neighbors *neighbors, Protein *protein);
int clean_short_structures(Ideal_struct * ideal_helix, Ideal_struct * ideal_strand, Protein * protein);
int calculate_neighbors_distances(Protein *protein, Neighbors * neighbors);
int structure2sse (Protein *protein);
void return_missing_data_distances(Neighbors * neighbors);
double calc_dist_res(Residue *res1, Residue * res2);
char struct_type(Residue * res);
int trace_back(int count, int count_min, int id, Protein * protein);
int is_regular_struct(Neighbors * neighbors, Ideal_struct * ideal_struct);
int enumerate_structures(Protein *protein);
int fit_line (double **point, int no_points, double center[3], double direction[3]);
int process_sse (Protein * protein, int type, double ** point, int number_of_points,
		 int first_res_index, int last_res_index, SSElement * element);
