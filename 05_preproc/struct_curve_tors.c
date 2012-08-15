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
 * Function that calculates the first determination f' of the function f
 * f = x[0] + x[1]*t + x[2]*t^2 + x[3]*t^3
 * f' = x[1] + 2 * x[2]*t + 3 * x[3]*t^2
 * @param x 
 * @param t
 * @return value of f'
 */

double der1(double *x, double t) {
    return x[1] + 2 * x[2] * t + 3 * x[3] * pow(t,2);
}


/**
 * Function that calculates the second determination f'' of the function f
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
 * Function that calculates the third determination f''' of the function f
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
int find_beta_curvature(Protein * protein){
    int i, j;
    int length = protein->length;
    double x[NO_OF_POINTS], y[NO_OF_POINTS], z[NO_OF_POINTS]; // coordinates of consecutive CA atoms
    Residue res;
    int exists_ca = 1;
    float curv;
    // printf("%d\n", length);
    
    for (i=0; i < length -(NO_OF_POINTS -1); ++i) {
       // printf("%d\n", i);
        exists_ca = 1;
        for (j=0; j<NO_OF_POINTS; ++j){
            res = protein->sequence[i+j];
            if ( ! res.Ca) {
                exists_ca = 0;
                break;
            }
            
            // taking coordinates
            x[j] = res.Ca->x;
            y[j] = res.Ca->y;
            z[j] = res.Ca->z;
        }
        if (! exists_ca) continue;
        // calculation of the curvature and the structure of the middle residue 
        curv = fit_curve(x,y,z); 
        if (curv < MAX_CURVE) {
            protein->sequence[i+NO_OF_POINTS/2].belongs_to_strand = 1;
            protein->sequence[i+NO_OF_POINTS/2].belongs_to_helix = 0;
        }
       // printf("%d %c %lf\n", i+1+2, struct_type(&protein->sequence[i+2]), curv);
    }
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

void curvature(double *a, double *b, double *c, double *t, int size, double *curv){
    double *num = (double *) malloc(size * sizeof(t));
    double *denum = (double *) malloc(size * sizeof(t));
    
    int i;
    for (i = 0; i < size; ++i) {
        num[i] = sqrt(pow(expr2(b,c,t[i]),2) + pow(expr2(c,a,t[i]),2) + pow(expr2(a,b,t[i]),2));
        denum[i] = pow((pow(der1(a, t[i]),2) + pow(der1(b, t[i]),2) + pow(der1(c, t[i]),2)), 1.5);
        curv[i] = num[i]/denum[i];
    }
    free(num);
    free(denum);
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

void torsion(double *a, double *b, double *c, double *t, int size, double *tors){
    double *num = (double *) malloc(size * sizeof(t));
    double *denum = (double *) malloc(size * sizeof(t));
    
    int i;
    for (i = 0; i < size; ++i) {
        num[i] = expr1(a,b,c,t[i]) +  expr1(b,c,a,t[i]) + expr1(c,a,b,t[i]);
        denum[i] = pow(expr2(b,c,t[i]),2) + pow(expr2(c,a,t[i]),2) + pow(expr2(a,b,t[i]),2);
        tors[i] = num[i]/denum[i];
    }
    free(num);
    free(denum);
}

/** Function that return init parameters for fitting. The initial guess is a line between the 
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
    
/*
    
    p0[2] = 0;
    p0[3] = 0;
    p0[6] = 0;
    p0[7] = 0;
    p0[10] = 0;
    p0[11] = 0;
*/
}


/**
 * function that we use for estimation of a curve that pass through five consecutive C alpha atoms
 * @param t
 * @param p
 * @return 
 */

double f( double t, const double *p )
{
    return p[0] + p[1]*t + p[2]*pow(t,2) + p[3]*pow(t,3);
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

double fit_curve(double *x, double *y, double *z){
    
    /* parameter vector */

    int n_par = 12; // number of parameters in model function f
    double par[12] = {0}; // arbitrary starting value
    
    /* data points */
    
    int m_dat = 3 * NO_OF_POINTS;

    /* auxiliary parameters */

    lm_status_struct status;
    lm_control_struct control = lm_control_double;
    //control.printflags = 3; // monitor status (+1) and parameters (+2)
    
    //control.maxcall = 10;
    //control.ftol = 0.00001;
    
    int i;
    
    double t[NO_OF_POINTS];
    
    double step = 1./(NO_OF_POINTS -1);
    
    for (i = 0; i < NO_OF_POINTS ; i++) {
        t[i] = i * step; 
    }
    
    srand((unsigned int)time(NULL));

    
/*
    for (i = 0; i < 5; ++i) {
        printf("%lf ", x[i]);
    }
    
    for (i = 0; i < 5; ++i) {
        printf("%lf ", y[i]);
    }
    
    for (i = 0; i < 5; ++i) {
        printf("%lf ", z[i]);
    }
*/
    
    set_init_param(x, y, z, par);
    data_struct data = { t, x, y, z, f };
    
/*
    for (i = 0; i < 12; ++i) {
        printf("%lf ", par[i]);
    }
    printf("\n");
*/
    
    /* perform the fit */

    //printf( "Fitting:\n" );
    
    
    
   
    lmmin( n_par, par, m_dat, (const void*) &data,
           residuals, &control, &status, lm_printout_std );

    double a[4], b[4], c[4];
    
    for (i = 0; i < 4; i++) {
        a[i] = par[i];
        b[i] = par[i+4];
        c[i] = par[i+8];
    }
    
     
     if (status.fnorm > 2) {
         return 1;
     }
  
/*
    for (i = 0; i < 12; ++i) {
        printf("%lf ", par[i]);
    }
    
    printf("\n");
*/
    
    double curv[NO_OF_POINTS];
    curvature(a, b, c, t, NO_OF_POINTS, curv);
    
    double output;
    output = curv[NO_OF_POINTS/2];
    
    return output; // return the curvature of the middle element
}
