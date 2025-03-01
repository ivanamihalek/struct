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


double int_table [TABLE_SIZE][TABLE_SIZE];

double exp_table [TABLE_SIZE];


double lookup ( double cos_alpha, double cos_beta ) {

    int a, b;

    a = (int) ( (1.0 - cos_alpha)*NR_POINTS/2 );
    b = (int) ( (1.0 - cos_beta)*NR_POINTS/2 );

    if (a > NR_POINTS) a = NR_POINTS; 
    if (a < 0) a = 0; 
    if (b > NR_POINTS) b = NR_POINTS; 
    if (b < 0) b = 0;
    
    return int_table[a][b];
}

int set_up_exp_table () {
    int i;
    for (i=0; i< TABLE_SIZE; i++)
	exp_table[i] = exp (-(double)i*MAX_EXP_VALUE/(TABLE_SIZE-1));

    return 0;
}
