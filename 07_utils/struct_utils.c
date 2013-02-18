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


void * emalloc(int  size)
{
    void * ptr;
    if ((ptr = calloc(size, 1)) == NULL) {
	fprintf (stderr,  "emalloc: no memory for %u bytes", size);
	return NULL;
    }

    return ptr;
}



FILE * efopen(char * name, char * mode) {

    FILE * fp;

    if ((fp = fopen(name, mode)) == NULL) {
	fprintf (stderr,  
	      "Cannot open \"%s\" for \"%s\"\n", name, mode);
	return NULL;
    }

    return fp;
}


/**********************************************************************/  
int infox ( char * errmsg, int exitval) {
    fprintf (stderr, "%s\nExiting at %s:%d.\n",		
	     errmsg, __FILE__, __LINE__ );		
    exit (exitval);
}


/**********************************************************************/  
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch]  */
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

float **fmatrix(int rows, int columns){
    float **m;
    int i;
        /* allocate pointers to rows */
    m=(float **) malloc(rows*sizeof(float*));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in chmatrix().\n");
	return NULL;
    } 
    /* allocate rows and set pointers to them */
    m[0]=(float *) calloc( rows*columns, sizeof(float));
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


double ***d3matrix(int rows, int columns, int floors){
    double ***m;
    int i, j;
    /* allocate pointers to rows */
    m=(double ***) malloc(rows*sizeof(double**));
    if (!m)  {
	fprintf (stderr,"row allocation failure  in d3matrix().\n");
	return NULL;
    }
    /* allocate rows and set pointers to them */
    m[0]=(double **) calloc( rows*columns, sizeof(double*));
    if (!m[0]) {
	fprintf (stderr,"column allocation failure in d3matrix().\n");
 	return NULL;
    }
    /* allocate overall block */
    m[0][0]=(double *) calloc( rows*columns*floors, sizeof(double));
    if (!m[0][0]) {
	fprintf (stderr,"column allocation failure in d3matrix().\n");
 	return NULL;
    }
    for( j=1; j < columns; j++) {
	m[0][j] = m[0][j-1] + floors;
    }
    for( i=1; i < rows; i++) {
	m[i]    = m[i-1] + columns;
	m[i][0] = m[i-1][0] + columns*floors;
	for( j=1; j < columns; j++) {
	    m[i][j] = m[i][j-1] + floors;
	}
    }
    /* return pointer to array of pointers to rows */ 
    return m; 
}



/* free a  matrix  */
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
void free_fmatrix(float **m)
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

int intmatrix_init(int **matrix, int rows, int columns, int val){
    int i,j;
    for (i=0; i<rows; i++) {
	for (j=0; j<columns; j++) {
	    matrix[i][j] = val;
	}
    }
    return 0;
}


/***********************************************************************/
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
/***************************************************************************/
int tokenize ( char token[MAX_TOK][MEDSTRING], int * max_token,
	       char * line , char comment_char) {
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
    while ( !string[ctr] && ctr < length) ctr++;
    
    if ( ctr == length ) return 1; /* empty string */
    
    if ( ctr ) {
	memmove (string, string+ctr, length-ctr);
	memset ( string+length-1-ctr, 0, ctr);
    }

    return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

int  cluster_counter (int  no_of_things,  int **neighbors, 
		      int * no_of_clusters, int ** clusters) {

    /* arrays */ 
    int*  * flag_ptr;          /* array of pointers to cluster flags*/
    int   * flags;             /* array of available flags */
    int   * bin  ;             /* for counting different clusters */
    /* counters, booleans etc */
    int flag_ctr, this_thing, other_thing; 
    int new_flag, cluster_count, isolated_count;
    int this_value, other_value;
    int color;
    
    /* flag ptrs   */
    flag_ptr     = emalloc (no_of_things*sizeof(int*));
    /* flags        */
    flags        = emalloc (no_of_things*sizeof(int));
    /* bins         */
    bin          = emalloc (no_of_things*sizeof(int));
    /* check if all alive: */ 
    if ( !( flag_ptr && flags && bin) ) return 1;
	
    /* set all the flags to 0 */
    memset (flags, 0, no_of_things*sizeof(int));
    /* the number of times new flag is assigned:*/
    new_flag = 0;
    /* set all the flag ptrs  to NULL */
    memset (flag_ptr, 0, no_of_things*sizeof(int*));
    /* color by cluster */ 
     
    for (this_thing=0; this_thing < no_of_things; this_thing++) {
	for (other_thing=this_thing+1; other_thing < no_of_things; other_thing++) {
	    if ( neighbors[this_thing][other_thing]) {
		if (flag_ptr[this_thing]){
		    if (flag_ptr[other_thing]){ /*if both ptrs assigned*/
			/*************************************************/
			/* look at the flag values they are assigned to: */
			if ( *flag_ptr[this_thing]  !=  *flag_ptr[other_thing] ) { 
			    /* i.e. do something only if they differ*/
			    this_value   = *flag_ptr[this_thing];
			    other_value  = *flag_ptr[other_thing];
			    for ( flag_ctr=0; flag_ctr < new_flag; flag_ctr++ ) {
				if ( flags[flag_ctr] == other_value) {
				    flags[flag_ctr] = this_value;
				}
			    }
				    
			}
		    } else {                       /* one not assigned*/ 
			/*************************************************/
			flag_ptr[other_thing] = flag_ptr[this_thing];
		    }
		} else {
		    if (flag_ptr[other_thing]){ /* one not assigned*/
			/*************************************************/
			flag_ptr[this_thing]  = flag_ptr[other_thing];
		    } else {                      /* both null*/
			/*************************************************/
			/*  create new flag*/
			flags[new_flag] = new_flag;
			/*  make both ptrs point there*/
			flag_ptr[this_thing]  = flag_ptr[other_thing] = &flags[new_flag];
			new_flag++;
		    }
		    
		}

	    }
	}
    }

    /*count the clusters*/
    memset (bin, 0, no_of_things*sizeof(int));
    for (cluster_count=0; cluster_count <=no_of_things; cluster_count++ ) {
	 memset (clusters[cluster_count], 0, (no_of_things+1)*sizeof(int));
    }
    cluster_count = 0;
    isolated_count = 0;
    for (this_thing=0; this_thing < no_of_things; this_thing++) {
	if ( !flag_ptr[this_thing] ) {
	    isolated_count++;	/* "cluster" 0 are isolated elements !*/
	    clusters [0][0]++;
	    clusters [0][ clusters [0][0] ] = this_thing;
	} else {
	    color = *flag_ptr[this_thing];
	    if ( ! bin[color] ){
		cluster_count ++;
	    }
	    bin[color] ++;
	    color += 1;
	    clusters [color][0]++;
	    clusters [color][ clusters [color][0] ] = this_thing;
	}
    }



    * no_of_clusters = cluster_count+isolated_count;
    
    free (flag_ptr); 
    free (flags);
    free (bin);   
    return 0;
 }



/*******************************************************************/
int improvize_name ( char *filename, char chain, char *outstring) {

    int name_length = strlen(filename);
    int c, tokenctr;
    char token[MAX_TOK][MEDSTRING];
    char scratchstring[MEDSTRING] = {'\0'};;
    char auxstring[MEDSTRING] = {'\0'};;
    int maxtoken;
    char comment_char = '!';

    /* this assumes that the scratchstring is empty; we'll not take care of that here */
    
    for (c=name_length-1; c>=0; c--) {
	if (filename[c] == '/') filename[c] = ' ';	
    }
    
    tokenize ( token, &maxtoken, filename, comment_char);

    sprintf ( auxstring, "%s", token[maxtoken]);

    for (c=strlen(auxstring)-1; c>=0; c--) {
	if (auxstring[c] == '.') auxstring[c] = ' ';	
    }
    tokenize ( token, &maxtoken, auxstring, comment_char);

    if (!strcmp(token[maxtoken], "pdb") || !strcmp(token[maxtoken], "ent") ) maxtoken --;

    sprintf (scratchstring, "%s", token[0]);
    for (tokenctr=1; tokenctr<=maxtoken; tokenctr++) {
	sprintf (scratchstring, "%s.%s", scratchstring, token[tokenctr]);
    }
    
    /*  */
    if ( !strncmp(scratchstring, "pdb",3) ) {
      memset (&auxstring[0], 0, MEDSTRING*sizeof(char));
      memcpy (&auxstring[0], &scratchstring[3], (strlen(scratchstring)-3)*MEDSTRING);
    }
    if (chain) scratchstring[strlen(scratchstring)] = chain;

    sprintf (outstring, "%s", scratchstring);
    
    return 0;
}
