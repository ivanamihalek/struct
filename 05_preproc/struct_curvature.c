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

#include "struct_lmmin.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"


/**
 * Function that calculates the first derivative  of the polynomial of the third order
 * f = x[0] + x[1]*t + x[2]*t^2 + x[3]*t^3
 * f' = x[1] + 2 * x[2]*t + 3 * x[3]*t^2
 * @param x 
 * @param t
 * @return value of f'
 */

double der1(double *x, double t) {
    return x[1] + 2 * x[2] * t + 3 * x[3] * t*t;
}


/**
 * Function that calculates the second derivative of the polynomial of the third order
 * f = x[0] + x[1]*t + x[2]*t^2 + x[3]*t^3
 * f'' = 2 * x[2] + 6 * x[3]*t
 * @param x
 * @param t
 * @return 
 */
double der2(double *x, double t) {
    return 2 * x[2] + 6 * x[3] * t;
}

/**
 * Function that calculates the third derivative f''' of the polynomial of the third order
 * f = x[0] + x[1]*t + x[2]*t^2 + x[3]*t^3
 * f''' = 6 * x[3]
 * @param x
 * @param t
 * @return 
 */

double der3(double *x, double t) {
    return 6 * x[3];
}

/**
 * Function that calculates an expression necessary for curvature and torsion calculation
 * @param a
 * @param b
 * @param t
 * @return 
 */

double expr2(double * a, double * b, double t){
    return der1(a,t)*der2(b,t) - der2(a,t)*der1(b,t);
}

/**
 * Function that calculates an expression necessary for curvature and torsion calculation
 * @param a
 * @param b
 * @param c
 * @param t
 * @return 
 */

double expr1(double * a, double * b, double *c, double t) {
    return der3(a, t)*expr2(b,c,t);
}


/**
 * Function that determines if the specific residue is the element of a beta strand.
 * The determination is based on curvature calculation
 * 
 * @param protein
 * @return 
 */
int strand_by_curvature(Protein * protein){
    int i, j;
    int length = protein->length;
    int no_nbrs = FITTING_NEIGHBORHOOD;
    int no_of_points = 2*FITTING_NEIGHBORHOOD +1;
    double x[no_of_points], y[no_of_points], z[no_of_points]; // coordinates of consecutive CA atoms
    Residue res;
    int exists_ca = 1, helix_nbr = 0;
    double a[4], b[4], c[4];

    // determining that something is a helix is easier. so we take this as done already
    // we do not mess (or waste time) with regions already established to be helices
    
    for (i=no_nbrs; i < length-no_nbrs; ++i) {

	//if (protein->sequence[i].belongs_to_helix) continue;

	exists_ca = 1;
	helix_nbr = 0;
        for (j=0; j< no_of_points; ++j){
            res = protein->sequence[i-no_nbrs+j];
            if ( ! res.Ca) {
                exists_ca = 0;
                break;
            }
            //if (protein->sequence[i-no_nbrs+j].belongs_to_helix)  {
	    //	helix_nbr = 1;
	    //	break;
	    //}
            //  coordinates
            x[j] = res.Ca->x;
            y[j] = res.Ca->y;
            z[j] = res.Ca->z;
        }
        if (! exists_ca || helix_nbr) continue;
        // use polynomial fit to estimate the curvature
        polynomial_fit(x,y,z, no_of_points, a, b, c);
	
	//printf("%3d   %8.3lf  %8.3lf  \n", i, curvature(a, b, c, 0.5), fabs(torsion(a, b, c, 0.5)));
 	
        //if (   curvature(a, b, c, 0.5) < MAX_CURVATURE && fabs(torsion(a, b, c, 0.5)) < MAX_TORSION) {
	if (   curvature(a, b, c, 0.5) < MAX_CURVATURE ) {
            protein->sequence[i].belongs_to_strand = 1;
            protein->sequence[i].belongs_to_helix  = 0;
        }
	
    }

    // second pass: check one more time whether the strand we are defining is too long:
    /* int candidate_strand_length = 0; */
    /* int strand_start = -1; */
    /* for (i=no_of_points/2; i < length -(no_of_points/2 -1); ++i) { */
    /* 	if (! protein->sequence[i].belongs_to_strand) { */

    /* 	    if ( candidate_strand_length >  MIN_CaS_IN_SSE) { */
    /* 		printf (" strand_start %3d  end %d \n",  strand_start+1, i); */
    /* 	    } */
	    
    /* 	    candidate_strand_length = 0; */
    /* 	    strand_start = -1; */
    /* 	    continue; */
    /* 	} */
    /* 	if ( strand_start < 0 ) strand_start = i; */
    /* 	candidate_strand_length ++; */
    /* } */
    return 0;
}


/**
 * Function that calculates the curvature of a curve using the following formula
 * 
 * K=(((z''y'-y''z')^2 + (x''z'-z''x')^2 + (y''x'-x''y')^2)^0.5)/(x'^2+y'^2+z'^2)^1.5
 * x(a,t) = a[0] + a[1]*t +a[2]*t^2 + a[3]*t^3
 * y(b,t) = b[0] + b[1]*t +b[2]*t^2 + b[3]*t^3
 * z(b,t) = c[0] + c[1]*t +c[2]*t^2 + c[3]*t^3
 * 
 * @param a parameter
 * @param b parameter
 * @param c parameter
 * @param t parameter
 * @param size number of points
 * @param curv output array of curvature values
 */

double curvature(double *a, double *b, double *c, double t){
    double num;
    double denum;
    double curv;
    num     = sqrt(pow(expr2(b,c,t),2) + pow(expr2(c,a,t),2) + pow(expr2(a,b,t),2));
    denum   = pow((pow(der1(a, t),2)   + pow(der1(b, t),2) + pow(der1(c, t),2)), 1.5);
    if ( denum < 1.0-10) return MAX_CURVATURE+1;
    curv    = num/denum;

    //printf ("   %8.3lf  %8.3lf     %8.3lf\n",  num, denum, curv); 

    return curv;
}

/**
 * Function that calculates the torsion of a curve using the following formula
 * T = (x'''(y'z''-y''z')+y'''(x''z-x'z'')+z'''(x'y''-x''y'))/((y'z''-y''z')+(x''z-x'z'')+(x'y''-x''y'))
 * 
 * x(a,t) = a[0] + a[1]*t +a[2]*t^2 + a[3]*t^3
 * y(b,t) = b[0] + b[1]*t +b[2]*t^2 + b[3]*t^3
 * z(b,t) = c[0] + c[1]*t +c[2]*t^2 + c[3]*t^3
 * 
 * @param a
 * @param b
 * @param c
 * @param t
 * @param size of the list of parameters t
 * @param tors
 */

double torsion(double *a, double *b, double *c, double t){
    double num;
    double denum;
    double tors;
    
    num   = expr1(a,b,c,t) +  expr1(b,c,a,t) + expr1(c,a,b,t);
    denum = pow(expr2(b,c,t),2) + pow(expr2(c,a,t),2) + pow(expr2(a,b,t),2);
    if ( denum < 1.0-10) return MAX_TORSION+1;
    tors  = num/denum;

    return tors;
    
}

/** Function that returns init parameters for fitting. The initial guess is a line between the 
 * first and the last point
 * 
 * @param x
 * @param y
 * @param z
 * @param p0 return value via a pointer
 */
void set_init_param(double *x, double *y, double *z, double *p0) {
    p0[0] = x[0];
    p0[1] = x[4] - x[0];
    p0[4] = y[0];
    p0[5] = y[4] - y[0];
    p0[8] = z[0];
    p0[9] = z[4] - z[0];
    
    
    p0[2] = 10*(rand()/(RAND_MAX + 1.0) - 0.5);
    p0[3] = 10*(rand()/(RAND_MAX + 1.0) - 0.5);
    p0[6] = 10 *(rand()/(RAND_MAX + 1.0) - 0.5);
    p0[7] = 10 *(rand()/(RAND_MAX + 1.0) - 0.5);
    p0[10] = 10 *(rand()/(RAND_MAX + 1.0) - 0.5);
    p0[11] = 10 *(rand()/(RAND_MAX + 1.0) - 0.5);
}


/**
 * function that we use for estimation of a curve that pass through five consecutive C alpha atoms
 * the so called polynomial, Mile
 * @param t
 * @param p
 * @return 
 */

double f( double t, const double *p )
{
    double tsq = t*t;
    return p[0] + p[1]*t + p[2]*tsq + p[3]*tsq*t;
}


/**
 * function evaluation, determination of residues 
 * @param par
 * @param m_dat
 * @param data
 * @param fvec
 * @param info
 */

void residuals( const double *par, int m_dat, const void *data,
                       double *fvec, int *info )

//void residuals( double *par, double *fvec, int m_dat, int n_dat, void *data)
{
    /* for readability, perform explicit type conversion */
    data_struct *mydata;
    mydata = (data_struct*)data;
   
    
    int i;
    double a[4], b[4], c[4];
    
    for (i = 0; i < 4; i++) {
        a[i] = par[i];
        b[i] = par[i+4];
        c[i] = par[i+8];
    }
    
    for ( i = 0; i < m_dat/3; i++ ) {
	fvec[i*3+0] = mydata->x[i] - mydata->f( mydata->t[i], a );
        fvec[i*3+1] = mydata->y[i] - mydata->f( mydata->t[i], b );
        fvec[i*3+2] = mydata->z[i] - mydata->f( mydata->t[i], c );
    }    
            
}

/**
 * Function that fit a 3D curve using Levenberg Marquardt Least Squares Fitting
 * and return curvature of the curve
 * Input are coordinates of 5 consecutive C alpha atoms
 * 
 * x(a,t) = a[0] + a[1]*t +a[2]*t^2 + a[3]*t^3
 * y(b,t) = b[0] + b[1]*t +b[2]*t^2 + b[3]*t^3
 * z(b,t) = c[0] + c[1]*t +c[2]*t^2 + c[3]*t^3
 * 
 * @param x
 * @param y
 * @param z
 * @return 
 */

// warpper for the function below - different format for the points
int polynomial_fit_pointlist(double **pointlist, int no_of_points, double a[], double b[], double c[]){
    double *x, *y, *z;
    x=emalloc(no_of_points*sizeof(double));
    y=emalloc(no_of_points*sizeof(double));
    z=emalloc(no_of_points*sizeof(double));
    int i;
    for (i=0; i<no_of_points; i++) {
	x[i] = pointlist[i][0];
	y[i] = pointlist[i][1];
	z[i] = pointlist[i][2];
    }
    polynomial_fit (x, y, z,  no_of_points, a, b, c);
    free(x);
    free(y);
    free(z);
    return 0;
}

int polynomial_fit(double *x, double *y, double *z, int no_of_points, double a[], double b[], double c[]){
    
    /* parameter vector */
    
    int n_par = 12; // number of parameters in model function f
    double par[12] = {0}; // arbitrary starting value
    
    /* data points */
    
    int m_dat = 3 * no_of_points;

    /* auxiliary parameters */

    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    //control.printflags = 3; // monitor status (+1) and parameters (+2)
    
    //control.maxcall = 10;
    //control.ftol = 0.00001;
    
    int i;
    
    double t[no_of_points];
    
    double step = 1./(no_of_points -1);
    
    for (i = 0; i < no_of_points ; i++) {
        t[i] = i * step; 
    }
    
    srand((unsigned int)time(NULL));//?

    set_init_param(x, y, z, par);
    data_struct data = { t, x, y, z, f };
    
   
    lmmin( n_par, par, m_dat, (const void*) &data,
           residuals, &control, &status, lm_printout_std );

    for (i = 0; i < 4; i++) {
        a[i] = par[i];
        b[i] = par[i+4];
        c[i] = par[i+8];
    }
    
    if (status.fnorm > 2) {//?
         return  1.0;
     }
     
     return  0; 
}
