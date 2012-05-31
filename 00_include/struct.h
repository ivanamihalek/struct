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

# ifndef _BC_H
# define _BC_H
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>
# include <sys/stat.h>
# include "struct_hungarian.h"
# include "struct_utils.h"
# include "struct_pdb.h"
# include "struct_geometry.h"
# include "struct_structure2sse.h"
# ifdef DMALLOC
#   include "dmalloc.h"
# endif


# define  BUFFLEN 150


/* input type */
# define PDB  1
# define  DB  2

# define INIT_ALLOC_N  120

/* used in preprocessor */
# define MIN_HELIX_LENGTH 8
# define MIN_STRAND_LENGTH 3

/* pdb parser errors */
# define ERR_SSE_NONE    2
# define ERR_CA_TRACE    4
# define ERR_NONSENSE    8
# define ERR_CONTAINER  16
# define ERR_OVERLAP    32
# define ERR_MAX_ATOMS  64
# define ERR_NO_FILE_OR_CHAIN  128

# define MAP_SIM_THRESHOLD 0.5

/* used in alignment functions */
# define FAR_FAR_AWAY -1000

/* for matching elements of different type */
#define REWARD   1
#define PENALTY  0

typedef enum   /* Declares an enumeration data type called ALGORITHM */
{
    SEQUENTIAL,     
    OUT_OF_ORDER,
    BOTH
} ALGORITHM; 


typedef struct {
    char * tgt_filename;  /* input filename for the target (set); pointer to argv, if given */
    char * qry_filename;  /* input filename for the qry (set) */

    char * tgt_db;        /* for testing purposes only - outside db, separate from the pdb */
    char * qry_db; 
    double merge_cosine;  /* min cos angle between two SSEs to be represented
			     by (merged into) the same vector */
    double alpha;         /* gaussian width for the scoring fn F */
                          /* note that it is set from the input table */
                          /* rather than from the cmd file        */
    double tol;           /* ? */
    double z_max_store;   /* z-score cutoff to store a map */
    double z_max_out;     /* z-score cutoff for map output  */
    double F_guess_max;   /* max value of F for an intial guess    */
    double F_eff_max;     /* max value of F after gradient descent */
    double z_max_corr;    /* max value of z for which two maps are */
			  /* considered correlated */
    double z_min_compl;   /* min value of z for which two maps are */
			  /* considered complementary */
    double grad_step_size;/* step size in gradient descent */
    double grad_stop_tol; /* stopping precision in gradient descent */
    double far_far_away;  /* nonsense distance for the pairwise alignment*/
    double gap_open;      /* gap penalties for the pairwise alignment */
    double gap_extend;
    double endgap;
    double threshold_distance; /* the max distance for exhaustive search
				  for init triple of SSEs */

    double distance_tol_in_bb_almt; /* exp fallof for the bb almt score */
    
    double far_away_cosine;       /* minimum cosine for F_effective estimate */
    double H_length_mismatch_tol; /* tolerance for the difference in length
				     in matched helices */
    double S_length_mismatch_tol; /* tolerance for the difference in length
				     in matched strands */
    double avg_length_mismatch_tol;
    
    int min_no_SSEs;      /* min number of SSEs we are willing to consider as a hit */
    int exhaustive;       /* try all SSE triples (within the threshold (above)
			     instead of consecutive only            */
    int smith_waterman;   /* use Smith-Waterman rather than Needleman-Wunsch        */
    int omp;              /* use omp parallelization of the exhaustive search */

    ALGORITHM search_algorithm;  /* database search algorithm  {out_of_order, sequential, both*/
    ALGORITHM current_algorithm; /* current database search algorithm {out_of_order, sequential} */
    
    int use_endgap;
    int grid_size;        /* minimal number of points for the sphere grid */
    int number_maps_cpl;  /* number of top scoring maps among which to   */
                          /* look for a complement */ 
    int number_maps_out;  /* number of top scoring maps to output     */
    int grad_max_step;    /* max number of steps in  gradient descent */
    int exp_table_size;   /* size of the lookup table for the exp funcion */
    int verbose;          /* toggles verbose output */
    int use_length;       /* use length of SSEs in the structure matching */
    int print_header;     /* print header in the short output file  */
    int report_no_sse_overlap; /*produce short output line even when
				 there is no overlap in the SSE type */
    int report_no_match;
    int postprocess;      /* produce an output for postprocessing */
    int optimize;         /* optimize backbone alignment */
    int preproc_only;
    char chain_tgt;
    char chain_qry;
  
    char outdir [BUFFLEN]; /* output directory - will not be crated if not present   */
    char outname[BUFFLEN];/* root name for the output file(s)   */
    char path   [BUFFLEN];   /* path to the integral table */
    
} Options;

typedef struct {
    double theta;
    double phi;
    void * subgrid;
} Point;




typedef struct {
    char name[SHORTSTRING];
    int no_of_residues, no_of_elements; 
    int no_of_helices, no_of_strands;
    SSElement  * element; /*list of secondary structure elements*/
} Descr;


typedef struct{
    double ** full;
    double ** cm; 
    double ** translation; /* this one will change
			      with changing CM of the subset */
    double *  transl_norm; /* this one will change
			      with changing CM of the subset */
    double origin[3];      /* origin point that traslation vectors refer to */
    int * full_type;
    int * length;
    int N_full;
    int full_no_of_strands;
    int *full_sheet_id;
    
    double ** compact;
    int *compact_type;
    int N_compact;
    int **represents;
    int * is_rep_by;
    
} Representation;


# define MAX_MULT_MATCH 100

# define MAP_MAX   100 /*number of maps to keep */

typedef struct {
    
    ////////////////////
    // general
    int size;                // max number of mapped elements (SSEs or Ca, depending where we use the structure
    int matches;             // the number ofactually mapped SSEs
    double q[4]; /* rotation (quaternion)  -- the rotation which results in this particular map */
    double T[3]; /* translation (once we believe we know it)  -- the translation (for the Ca level)*/
    
    ////////////////////
    // referring to SSE-level map:
    int *x2y, *y2x;          // x2y: for each mapped pair (x,y), returns y, the index of the SSE in the structure Y,
                             // given x, the index of SSE in structure X  ( x2y[x] = y, e.g. x2y[3] = 7)
    int x2y_size, y2x_size;  // the size of the two maps above; ultimately the two should be equal
                             // in some general intermediate step that might not be the case, in the
                             // current implementation it always is
    double avg_length_mismatch; // average difference in the length of mapped SSEs
    double rmsd;             /* rmsd for the centers of matched SSEs */

    ////////////////////
    // "urchin" scoring
    double **cosine;          // table of angle cosines for all pairs (x,y) (all of them, not just the mapped ones)
    double **sse_pair_score;           // table of exp terms for all pairs (x,y) (all of them, not just the mapped ones)
    double F;                 // value of the  function F for this map
    double avg, avg_sq;       // refers to the avg and average square of F over the sapce of all rotations
                              // the input for the calulcation of the z-score
    double z_score;           // z-score for F (based on avg and avg_sq0
    double assigned_score;    // sum of the exp terms of but only for the  matched SSE pairs (x,y)
    
    ////////////////////
    // referring to Ca-level map:
    int *x2y_residue_level, *y2x_residue_level;  // the same as the x2y above, this time not on SSE, but on Ca level
    int x2y_residue_l_size, y2x_residue_l_size;
    int res_almt_length;    // length of the alignment on the Ca level
    double res_rmsd;        /* rmsd for the matched Ca atoms*/
    double aln_score;         // like the "assigned score" above, but for the Ca level
    double res_almt_score;/*not sure - cou;ld be I am duplicating aln_score */
    
    ////////////////////
    // complementary or sub-maps - never mind for now, just leave as is
    int *submatch_best;         // the best map which complements this one; don't worry about it right now
    double score_with_children; // this goes with the submatch above = never mind
    double compl_z_score;       // z-score for the submatch

    ///////////////////
    // file to which the corresponding pdb was written
    char filename[MEDSTRING];
    
} Map;


typedef struct {

    Map * map;
    int no_maps_allocated;
    int no_maps_used;
    
    int *map_best;
    int best_array_allocated;
    int best_array_used;
    
    int NX_allocated;
    int NY_allocated;

} List_of_maps;


typedef struct {
    double total_assigned_score;
    double fraction_assigned;
    double z_score;
    double w_submap_z_score;
    double gap_score;
    double rmsd;
    double conformer;
    double avg_length_mismatch;
    double submap_avg_length_mismatch;
    /* scores that appear in potprocessing */
    double res_almt_score;
    double res_rmsd; /* rmsd for aligned CAs */
    int res_almt_length;
    int number_of_maps;
    double q_rmsd[3];
} Score;

/*************************************/
/* penalty scheme for the alignment  */
/*************************************/
typedef struct {
    int endgap_special_treatment;
    double gap_opening;
    double gap_extension;
    double endgap;
    double *custom_gap_penalty_x; /* special, per-element gap penalties*/
    double *custom_gap_penalty_y;

} Penalty_parametrization;



/******************************/
/* integral lookup table :    */
/******************************/

# define NR_POINTS 100
# define TABLE_SIZE NR_POINTS+1

# define MAX_EXP_VALUE 10

extern double int_table [TABLE_SIZE][TABLE_SIZE];
extern double exp_table [TABLE_SIZE];
extern Options options; /* defined in struct.c */

/******************************/
/* function declarations :    */
/******************************/

double alignment_score (Protein * protein1, Protein * protein2, int * residue_map_i2j,
			    double **R, double *T,  double d0 );
int align_backbone (Descr *tgt_descr, Protein *tgt_structure, Representation * tgt_rep,
		    Descr *qry_descr, Protein *qry_structure, Representation * qry_rep, List_of_maps *list);
int check_gap_lengths (Map * map, double *gap_score );
int check_input_type  (FILE *fptr);
int close_digest      (clock_t CPU_time_begin, clock_t CPU_time_end, FILE *digest);
int compare_descr     (Descr *descr1, Descr *descr2, List_of_maps *list, Score *score);
int complement_match  (Representation* X_rep, Representation* Y_rep, List_of_maps * list);
int construct_translation_vecs (Representation *X_rep,  Representation *Y_rep, Map *map);
int descr_init ( Descr * description);
int descr_out (FILE * fptr, Descr * descr);
int descr_shutdown ( Descr * description );

double F (double **X, int * x_type,  int NX,
	  double **Y, int * y_type, int NY, double alpha);
int F_moments (double **x, int *x_type, int NX,
	       double **y, int *y_type, int NY, double alpha,
	       double *avg_ptr, double *avg_sq_ptr,
	       double *std_ptr);
int fill_protein_info ( FILE * fptr,  char chain, Protein * protein);
int find_next_pair (double **X, double **Y, int *x_type, int *y_type,
		    int NX, int NY, int nbr_range, int x_ctr, int y_ctr,
		    int *x_next_ptr, int *y_next_ptr);
int find_quat (double ** Xin, double **Yin, int no_vectors,
		 double * quat, double *qof);
int find_quat_exp (double ** X, int NX, double **Y, int NY,
		     double alpha, double * quat, double *qof);
int find_uniq_maps (List_of_maps  *list1, List_of_maps *list2, List_of_maps  *list_uniq);
int free_map (Map *map);
int geometric_descriptors (Protein * protein);
int get_next_descr (int input_type, FILE * fptr,  char chain, Protein *protein, Descr * description);
int hungarian_alignment (int NX, int NY, double **similarity, int * x2y, int * y2x, double * alignment );
int init_digest (Descr *qry_descr, Descr *tgt_descr, FILE ** digest);
int initialize_map (Map *map, int NX, int NY );
int input  (FILE * fptr, Descr * description);
int list_alloc (List_of_maps * list, int NX, int NY, int fake);
int list_shutdown (List_of_maps * list, int fake);
double lookup ( double alpha, double beta);
int map_assigned_score (Representation *X_rep,  Map* map);
int map_complementarity (Map *map1, Map *map2, double *z);
int map_consistence (int NX, int NY, Map *combined_map, Map *map1, Map *map2,
		     double *total_ptr, double * gap_score);
int map_reduced_reps (Representation *rep1, Representation *rep2, List_of_maps *list);
int mat_out (double A[4][4], char *name);
int mat_mult (double new [4][4], double  A[4][4], double  B[4][4]);
int mat_sum (double sum[4][4], double  new_term[4][4]);
int mat_diag ( double B[4][4], double eval[4], double evect[4][4] );
int mat_exp (double expB[4][4], double B[4][4]);
int match_length (int N, int *x2y);
int multiply (double *quat1, double *quat2_in,
	      int conjugate,  double *product);
int needleman_wunsch (int max_i, int max_j, double **distance,
		      int *map_i2j, int * map_j2i, double *aln_score);
int normalized_cross (double *x, double *y, double * v, double *norm_ptr);
int optimize_backbone_alignment (Descr *descr1, Protein * protein1, Representation *rep1, 
				 Descr *descr2, Protein * protein2, Representation *rep2, 
				 List_of_maps *list);
int output (FILE * fptr,char *name, char chain,  Protein * protein);
int print_map (FILE *fptr, Map * map, Descr * descr1, Descr * descr2,
	       Protein *protein1, Protein *protein2,  int tab);
    
double qnorm (double *quat);
double q_dotprod (double *quat1, double *quat2);
int quat_to_R ( double quat[4], double **R);
double quat_rmsd (double parent_map_q[4], double  current_map_q[4]);
int random_q ( double exp_s[4],  double theta_step );
int read_integral_table (char * file_name );
int read_pdb (char * pdbname,  char chain, Protein * protein);
int rep_initialize (Representation * rep,  Descr * descr  );
int rep_shutdown (Representation * rep);
int results_out  (Descr *tgt_descr, Protein *tgt_structure, Representation * tgt_rep,
		  Descr *qry_descr, Protein *qry_structure, Representation * qry_rep,
		  List_of_maps *list, FILE * digest);
int rotate(double **Ynew, int NY, double **R, double ** Y);
int set_match_algebra ();
int set_up_exp_table ();
void similarity_to_scoring(int NX, int NY, int multiplier, double** similarity, int ** hungarian_alignment ) ;
char single_letter ( char code[]);
int sse2descriptor (Protein *protein, Descr *descr);
int store_sse_pair_score (Representation *X_rep,  Representation *Y_rep, 
		 double **R, double alpha,  Map *map);
int structure2sse  (Protein *protein);


int smith_waterman (Penalty_parametrization *params, int max_i, int max_j, double **similarity,
		    int *map_i2j, int * map_j2i, double * aln_score);

int unnorm_dot (double *x, double *y, double * dot);
int vec_out (double *vec, int dim,  char * name );
int write_alignment (Protein *protein1, Protein *protein2, List_of_maps * list, Descr *descr1, Descr *descr2);
int write_digest(Descr *qry_descr, Descr *tgt_descr,
		 Representation * qry_rep, Representation * tgt_rep,
		 List_of_maps *list,FILE * digest);
int write_maps (FILE * fptr, Descr *descr1, Descr *descr2, List_of_maps *list);
int write_tfmd_pdb ( Protein * tgt_protein, List_of_maps *list, Descr *tgt_descr, Descr *qry_descr);
int find_Calpha (Protein *protein, int  resctr, double ca[3] );
double two_point_distance (double point1[3], double point2[3]);
int point_rot_tr (double point_in[3], double **R, double T[3],double point_out[3]);

# endif
