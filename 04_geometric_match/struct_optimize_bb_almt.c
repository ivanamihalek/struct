/*
This source code is part of deconSTRUCT,
protein structure database search and backbone alignment application.
Written by Ivana Mihalek, with contributions from Mile Sikic.
Copyright (C) 2012 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or,
at your option, any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# include "struct.h"

/**************************************************************/
/*  optimize the alignment on the backbone level, for each    */
/*  rotation+translation given in the map list                */
/**************************************************************/
int single_map_optimize_bb_almt (Protein * protein1, Protein * protein2, Map * map);

int optimize_backbone_alignment (Descr *descr1, Protein * protein1, Representation *rep1, 
				 Descr *descr2, Protein * protein2, Representation *rep2, 
				 List_of_maps *list){
    
     if ( list->no_maps_used == 0) return 1;
     
     int map_ctr,retval;
     
     Map *current_map;
       
     for (map_ctr=0; map_ctr<list->no_maps_used; map_ctr++) {
	 current_map = list->map+map_ctr;
	 retval = single_map_optimize_bb_almt (protein1, protein2, current_map);
	 if (retval) {
	     printf (" error optimize bb alignment   db:%s  query:%s \n", descr1->name, descr2->name);
	     exit (retval);
	 }
     }

     return 0;
}


int single_map_optimize_bb_almt (Protein * protein1, Protein * protein2, Map * map) {


    /* the current best guess for the rotation and translation are in
       map->q (the rotation representeed as a 4-component quaternion; the components defined as doubles),
       and  map->T (3 component; double); to get the rotation matrix use 
       quat_to_R (q, R); defined in  04_geometric_match/struct_quaternion.c:30
       Map is defined in 00_include/struct.h:190
       
    */


    /* after optimization replace map-> and map->T with the new values */
 
    
    return 0;
}
