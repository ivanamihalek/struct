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

/*****************************************************************************/
/*****************************************************************************/
int map_reduced_reps (Representation *rep1, Representation *rep2, List_of_maps *list) {
    
    /* descr 1 is db or target, descr2 is query - important in postprocessing*/
    int retval;
    int map_ctr;
    int NX, NY;

    /******************************************/
    /*  initialization                        */
    /******************************************/
    /*  shorthands   */
    NX     = rep1->N_full;
    NY     = rep2->N_full;
 
   
    if ( !list  || ! list->map )  return 1;
   
    
    /* the size has increased case */ 
    if ( NX > list->NX_allocated || NY > list->NY_allocated  ) {
	for ( map_ctr= 0; map_ctr< list->no_maps_allocated; map_ctr++) {
	    if ( list->map) if ( free_map( list->map+map_ctr) ) return 1;
	    
	    int NXalloc = NX > list->NX_allocated ? NX : list->NX_allocated;
	    int NYalloc = NY > list->NY_allocated ? NY : list->NY_allocated;
	    
	    if ( initialize_map( list->map+map_ctr, NXalloc, NYalloc) ) return 1;
	}
	list->NX_allocated = NX;
	list->NY_allocated = NY;
	
    } else { /* I should perhaps clean up anyway */
	for ( map_ctr= 0; map_ctr< list->no_maps_used; map_ctr++) {
	    clear_map (list->map+map_ctr);
	}
	list->no_maps_used = 0;
    }

    for ( map_ctr= 0; map_ctr<list->no_maps_allocated; map_ctr++) {
	(list->map+map_ctr)->x2y_size = NX;
	(list->map+map_ctr)->y2x_size = NY;
    }
    

    /******************************************/
    /*  look for maps between rep1 and rep2,  */
    /*    and for their complements           */
    /******************************************/
    retval  = complement_match (rep1, rep2, list);
    if (retval) return retval;
    

   
    return 0;
    
}

