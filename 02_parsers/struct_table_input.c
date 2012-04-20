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

int read_integral_table (char * filename ) {

    int a, b;
    int line_ctr;
    double value, delta = 0.0;;
    char line[BUFFLEN];
    FILE * fptr;
    /* int_table defined in struct_table.c and declared in struct.h */
    memset (int_table[0], 0, TABLE_SIZE*TABLE_SIZE*sizeof(double) );

    fptr = efopen (filename, "r" ) ;
    if ( !fptr ) return 1;

    line_ctr = 0;
    memset (line,  0, BUFFLEN);
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	line_ctr ++;
	if ( line[0] == '%' ) { /* comment */
	    if (strstr (line, "delta" ) ) {
		sscanf (line, " %*s  %*s   %lf", &delta );
	    }
	    continue;
	}
	if ( sscanf (line, "%d %d  %*s   %*s   %lf ",
		     &a, &b, &value ) == 3 ) {
	    if ( a >= TABLE_SIZE ||  b >= TABLE_SIZE  ) {
		fprintf (stderr, "Error reading in the lookup table:\n");
		fprintf (stderr, "line %d: %s\n", line_ctr, line);
		return 1;
	    }
	    int_table [a][b] = value;
	    int_table [b][a] = value;
	}
	memset (line,  0, BUFFLEN);
    }
    fclose (fptr);

    if ( !delta ) {
	fprintf (stderr, "Error reading in the lookup table:\n");
	fprintf (stderr, "no delta found\n");
	return 1;
    }

    options.alpha = delta;
    
    return 0;
    
}
