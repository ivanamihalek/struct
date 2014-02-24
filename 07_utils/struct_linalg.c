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



/************************************/
/************************************/
/************************************/
/************************************/
int sum_to_zero (double *a, double *b); /* defined below */
int parallel    (double *a, double *b); /* defined below */

/************************************/
int  unnorm_cross (double *x, double *y, double * v) {

    /* v is the output */

    if (parallel (x, y) ) return 1;
    
    v[0] = x[1]*y[2] -  x[2]*y[1];
    v[1] = x[2]*y[0] -  x[0]*y[2]; 
    v[2] = x[0]*y[1] -  x[1]*y[0];

    return 0;
    
}
/************************************/
int  normalized_cross (double *x, double *y, double * v, double *norm_ptr) {

    /* v is the output */
    double norm = 0;
    int i;

    if (parallel (x, y) ) return 1;
    
    v[0] = x[1]*y[2] -  x[2]*y[1];
    norm += v[0]*v[0];
    v[1] = x[2]*y[0] -  x[0]*y[2]; 
    norm += v[1]*v[1];
    v[2] = x[0]*y[1] -  x[1]*y[0];
    norm += v[2]*v[2];
    norm = sqrt (norm);

    for (i=0; i<3; i++ ) {
	v[i] /= norm;
    }

    if ( norm_ptr) *norm_ptr = norm;

    return 0;
    
}

/************************************/
int unnorm_dot (double *x, double *y, double * dot) {


    double cosine = 0;
    int i;
    
    for (i=0; i<3; i++ ) {
	cosine += x[i]*y[i];
    }

    if (cosine > 1.0 ) 
	cosine = 1.0; /* this should be numerical */
    if (cosine < -1.0 ) 
	cosine = -1.0; /* this should be numerical */
   
    *dot = cosine;
   
    return 0;
    
}

/***************************************/
int rotate_single_vector (double **R, double *x, double *x_rotated) {
    int i, j;
    for (i=0; i<3; i++ ) {
	x_rotated[i] = 0;
	for (j=0; j<3; j++ ) {
	    x_rotated[i] += R[i][j]*x[j];
	}
    }
    return 0;
}

/***************************************/
int rotate(double **Ynew, int NY, double **R, double ** Y) {

    int n, i,j;

    for (n=0; n<NY; n++ ) {

	for (i=0; i<3; i++ ) {
	    Ynew[n][i] = 0;
	    for (j=0; j<3; j++ ) {
		Ynew[n][i] += R[i][j]*Y[n][j];
	    }
	}
    }

    return 0;
}

/***************************************/
/***************************************/
/***************************************/
/*****************************************/
double norm (double *vec, int dim) {
    int i;
    double norm = 0;
    for (i=0; i<dim; i++) {
	norm +=  vec[i]* vec[i];
    }
    norm = sqrt(norm);
    return norm;

}

/*****************************************/
int vec_out (double *vec, int dim,  char * name ) {
    int i;
    double norm = 0;
    
    printf (" %5s  ", name);
    for (i=0; i<dim; i++) {
	printf ("%6.3lf  ", vec[i]);
	norm +=  vec[i]* vec[i];
    }
    printf ("   norm: %7.4lf  ", norm);
    printf ("\n");
    return 0;
}


/************************************/
int  sum_to_zero (double *a, double *b ) {
    int i;
    double sum,aux;

    sum = 0;
    for (i=0; i<3; i++ ) {
	aux = a[i] + b[i];
	sum += aux*aux;
    }
    return (sum < 0.001);
}

/************************************/
int  parallel (double *a, double *b ) {
    int i;
    double sum, diff, aux;

    sum  = 0;
    diff = 0;
    for (i=0; i<3; i++ ) {
	aux = a[i] + b[i];
	sum += aux*aux;
	
	aux = a[i] - b[i];
	diff += aux*aux;
    }
    return (sum < 0.00001 || diff < 0.00001);
}
