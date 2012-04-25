# ifndef _UTILS_H
# define _UTILS_H
# include <stdio.h>

int      array_qsort (int * sorted_pos, double * sa, int sequence_length );
char   **chmatrix(int rows, int columns);
double **dmatrix(int rows, int columns);
void *   emalloc(int	size);
FILE *   efopen(char * name, char * mode);
void     free_matrix(void **m);
void     free_cmatrix(char **m);
void     free_imatrix(int **m);
void     free_dmatrix(double **m);
void     free_d3matrix(double ***m);
void     free_strmatrix(char ***m);
int    **intmatrix(int rows, int columns);
char     single_letter ( char code[]);
int      string_clean ( char* string, int length);

# define MAX_EXP_VALUE 10
# define NR_POINTS 500
# define TABLE_SIZE NR_POINTS+1
extern double exp_table [TABLE_SIZE];
int set_up_exp_table ();
# endif
