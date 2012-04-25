# include "nossa.h"
#include <errno.h>
#include <fcntl.h>

#include <sys/types.h>

double exp_table [TABLE_SIZE];
int set_up_exp_table () {
    int i;
    for (i=0; i< TABLE_SIZE; i++)
	exp_table[i] = exp (-(double)i*MAX_EXP_VALUE/(TABLE_SIZE-1));

    return 0;
}

void * emalloc(int	size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}



FILE * efopen(char * name, char * mode)
{

    FILE * fp;


    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }

    return fp;

}




/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] 
 */
char **chmatrix(int rows, int columns){
    char **m;
    int i;
        /* allocate pointers to rows */
    m=(char **) malloc(rows*sizeof(char*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(char *) calloc( rows*columns, sizeof(char));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

int **intmatrix(int rows, int columns){
    int **m;
    int i;
        /* allocate pointers to rows */
    m=(int **) malloc(rows*sizeof(int*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(int *) calloc( rows*columns, sizeof(int));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}

double **dmatrix(int rows, int columns){
    double **m;
    int i;
        /* allocate pointers to rows */
    m=(double **) malloc(rows*sizeof(double*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(double *) calloc( rows*columns, sizeof(double));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in chmatrix().\n");
 	return NULL;
    }
    for( i=1; i < rows; i++)  m[i] = m[i-1] + columns;
    /* return pointer to array of pointers to rows */ 
    return m; 
}



/* free a  matrix  */
void free_matrix(void **m)
{
    free(m[0]);
    free(m);
}


/* free a  matrix  - save some warnings*/
void free_cmatrix(char **m)
{
    free(m[0]);
    free(m);
}
void free_imatrix(int **m)
{
    free(m[0]);
    free(m);
}
void free_dmatrix(double **m)
{
    free(m[0]);
    free(m);
}
void free_d3matrix(double ***m)
{
    free(m[0][0]);
    free(m[0]);
    free(m);
}
void free_strmatrix(char ***m)
{
    free(m[0][0]);
    free(m[0]);
    free(m);
}



/********************************************************************************************/


/* sort array according to the score in the other */
/* I couldn't declare pos_cmp within array_qsort  bcs it   crashed on mac */

double * score_array;

int pos_cmp (const void * a0, const void * b0) {
    
    int * a= (int*) a0;
    int * b= (int*) b0;
    if ( score_array[*a] > score_array[*b]) {
	return 1;
    }
    if ( score_array[*a] < score_array[*b]) {
	return -1;
    }
    return 0;
}

int array_qsort (int * sorted_pos, double * sa, int sequence_length ) {
    /* position comparison function */
    score_array = sa;

    qsort (sorted_pos, sequence_length, sizeof(int), pos_cmp);

    return 0;
}
/********************************************************************************************/
int tokenize ( char token[MAX_TOK][MEDSTRING], int * max_token, char * line , char comment_char) {
    /* assumes the tokens to be no bigger than MEDSTRING */ 
    
    char * chrptr, *last; 
    int current_token, current_char = 0;
    int reading;
   
    memset (token[0], 0, MAX_TOK*MEDSTRING*sizeof(char)); 
    chrptr = line;
    last   = chrptr + strlen (line);
    current_token = -1;
    current_char  =  0;
    reading = 0;
    while ( chrptr <= last) {
	if ( *chrptr == comment_char ) break;
	if ( *chrptr == '\n' ) break;
	if ( *chrptr && ! isspace(*chrptr) ) {
	    if ( ! reading ) {
		reading = 1;
		current_char = 0;
		current_token++;
		if ( current_token >= MAX_TOK ) {
		    return TOK_TOOMNY; /* defined in possum_utils.h */
		}
	    }
	    if ( current_char >= MEDSTRING ) {
		return TOK_TOOLONG;
	    }
	    token[current_token][current_char] = *chrptr;
	    current_char++;
	} else {
	    if ( reading ) {
		reading = 0;
	    }
	}
	chrptr++;
    }
    *max_token = current_token;

    return 0;
    
}

/**********************************************************/
/* get rid of spaces in a string */
int  string_clean ( char* string, int length) {
    int ctr;
    for (ctr = 0; ctr < length; ctr ++) {
	if ( isspace (string[ctr]) ) string[ctr] = '\0';
    }
    ctr=0;
    while ( !string[ctr]) ctr++;
    if ( ctr ) {
	memmove (string, string+ctr, length-ctr);
	memset ( string+length-1-ctr, 0, ctr);
    }

    return 0;
}

/**********************************************************/
/* three letter to single letter aa code                  */
char single_letter ( char code[]){

    switch ( code[0] ) {
    case 'A':
	switch ( code [1]) {
	case 'L':
	    return 'A';
	    break;
	case 'R':
	    return 'R';
	    break;
	case 'S':
	    switch ( code[2] ) {
	    case 'N':
		return 'N';
		break;
	    case 'P':
		return  'D';
		break;
	    }
	    break;
	}
	break;
    case 'C':
	return 'C'; 
	break;
    case 'G':
	/* the second letter is always L */ 
	switch ( code[2] ) {
	case 'U':
	    return 'E';
	    break;
	case 'N':
	    return  'Q';
	    break;
	case 'Y':
	    return 'G';
	    break;
	}
	break;
    case 'H':
	return  'H';
	break;
    case 'I':
	return  'I';
	break;
    case 'L':
	switch ( code [1]) {
	case 'E':
	    return 'L';
	    break;
	case 'Y':
	    return 'K';
	    break;
	}
	break;
    case 'M':
	return 'M';
	break;
    case 'P':
	switch ( code [1]) {
	case 'H':
	    return 'F';
	    break;
	case 'R':
	    return 'P';
	    break;
	}
	break;
    case 'S':
	return 'S';
	break;
    case 'T':
	switch ( code [1]) {
	case 'H':
	    return 'T';
	    break;
	case 'R':
	    return 'W';
	    break;
	case 'Y':
	    return 'Y';
	    break;
	}
	break;
    case 'V':
	return 'V';
	break;
	
    }


    fprintf (stdout, "Unrecognized amino acid code: %s.\n", code);
    return 0;
}



