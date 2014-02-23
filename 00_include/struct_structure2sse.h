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

int determine_sec_structure(Neighbors *neighbors, Protein *protein, int beta_curvature);
int clean_short_structures(Ideal_struct * ideal_helix, Ideal_struct * ideal_strand, Protein * protein);
int calculate_neighbors_distances(Protein *protein, Neighbors * neighbors);
int structure2sse (Protein *protein);
void return_missing_data_distances(Neighbors * neighbors);
double calc_dist_res(Residue *res1, Residue * res2);
char struct_type(Residue * res);
int trace_back(int count, int count_min, int id, Protein * protein);
int is_regular_struct(Neighbors * neighbors, Ideal_struct * ideal_struct);
int enumerate_SSEs(Protein *protein, int * counter);
int fit_line (double **point, int no_points, double center[3], double direction[3]);
int process_sse (Protein * protein, int type, double ** point, int number_of_points,
		 int first_res_index, int last_res_index, SSElement * element);
