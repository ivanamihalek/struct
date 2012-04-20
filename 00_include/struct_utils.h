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

# ifndef _UTILS_H
# define _UTILS_H


/******************************/
/*   tokenizer                */
/******************************/

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25


int  cluster_counter (int  no_of_things, int **neighbors, 
		      int *no_of_clusters, int **clusters);
int tokenize ( char  token[MAX_TOK][MEDSTRING], int * max_token,
	       char * line, char comment_char);




int       array_qsort (int * sorted_pos, double * sa, int sequence_length );
char    **chmatrix(int rows, int columns);
double  **dmatrix(int rows, int columns);
double ***d3matrix(int rows, int columns, int floors);
void *   emalloc(int	size);
FILE *   efopen(char * name, char * mode);
void     free_matrix(void **m);
void     free_cmatrix(char **m);
void     free_imatrix(int **m);
void     free_dmatrix(double **m);
void     free_d3matrix(double  ***m);
int    **intmatrix(int rows, int columns);
int intmatrix_init(int **matrix, int rows, int columns, int val);
int improvize_name ( char *filename,  char chain, char *outstring);
int      string_clean ( char* string, int length);
# endif
