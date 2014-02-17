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

# define MAX_CURVATURE  0.12
# define MAX_TORSION    1.2
# define FITTING_NEIGHBORHOOD  2 // number of points on each side we are fitting on
	 // 2*FITTING_NEIGHBORHOOD +1; must be greater than 4 bc we are fitting to poly of third degree
double der1(double *x, double t);
double der2(double *x, double t);
double der3(double *x, double t);

double expr2(double * a, double * b, double t);
double expr1(double * a, double * b, double *c, double t);

double curvature(double *a, double *b, double *c, double t);
double torsion(double *a, double *b, double *c, double t);

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

int strand_by_curvature(Protein * protein);
int polynomial_fit_pointlist(double **pointlist, int no_of_points, double a[], double b[], double c[]);
int  polynomial_fit(double *x, double *y, double *z, int no_of_points, double a[], double b[], double c[]);
#endif	/* STRUCT_CURVE_TORS_H */

