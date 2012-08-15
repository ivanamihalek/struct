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

#ifndef STRUCT_CURVE_TORS_H
#define	STRUCT_CURVE_TORS_H

# define MAX_CURVE 0.12
# define NO_OF_POINTS 5 // number of input points. Have to be >= 4

double der1(double *x, double t);
double der2(double *x, double t);
double der3(double *x, double t);

double expr2(double * a, double * b, double t);
double expr1(double * a, double * b, double *c, double t);

void curvature(double *a, double *b, double *c, double *t, int size, double *curv);
void torsion(double *a, double *b, double *c, double *t, int size, double *tors);

void set_init_param(double *x, double *y, double *z, double *p0); 
  
double f( double t, const double *p );

typedef struct {
    double *t;
    double *x;
    double *y;
    double *z;
    double (*f)( double t, const double *p );
} data_struct;

/* function evaluation, determination of residues */

void residuals( const double *par, int m_dat, const void *data,
                       double *fvec, int *info );

//void residuals( double *par, double *fvec, int m_dat, int n_dat, void *data);

double fit_curve(double *x, double *y, double *z);
int find_beta_curvature(Protein * protein);

#endif	/* STRUCT_CURVE_TORS_H */

