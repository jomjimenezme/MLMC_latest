/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Gustavo Ramirez, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori, Tilmann Matthaei, Ke-Long Zhang.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#include "main.h"

// vector storage for PRECISION precision
void vector_PRECISION_define( vector_PRECISION phi, complex_PRECISION value, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = value;
  } else {
    error0("Error in \"vector_PRECISION_define\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_define_random( vector_PRECISION phi, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = (PRECISION)(((double)rand()/(double)RAND_MAX))-0.5 + ( (PRECISION)((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I;
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );
}

// TODO: implement separate 3D funciton, this is just for testing
void vector_PRECISION_define_random_rademacher( vector_PRECISION phi, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
  for (int i = start; i < end; i++) {
    phi[i] = 0.0;
  }

  if(thread == 0){
    if ( phi != NULL ) {
      int i_global, i_local, site_index, owner;
      int dof_per_site = l->num_lattice_site_var; // variables per site
      int time_slice = 3;    // TODO: set from .ini file
      
      int depth = l->depth;
      int global_sites = g.global_lattice[depth][0] *
      g.global_lattice[depth][1] *
      g.global_lattice[depth][2] *
      g.global_lattice[depth][3];
      
      //  degrees of freedom in the full global vector
      int total_dof = dof_per_site * global_sites;
      
      // Number of dofs owned by each process //TODO: is always even?
      int chunk = total_dof / g.num_processes;
      
     
      // Loop over spatial coordinates 
      for (int z = 0; z < g.global_lattice[0][1]; z++) {
        for (int y = 0; y < g.global_lattice[0][2]; y++) {
          for (int x = 0; x < g.global_lattice[0][3]; x++) {
            
            // Get the global site index at (t_slice, z, y, x)
            site_index = lex_index(time_slice, z, y, x, g.global_lattice[0]);
            
            // Loop over the dof of that site
            for (i_global = dof_per_site * site_index;
                i_global < dof_per_site * site_index + dof_per_site;
            i_global++) {
              
              // Compute owner of the entry and local index 
              owner = i_global / chunk;
              i_local = i_global % chunk;
              
              // Only the owner assigns the value
              if(owner == g.my_rank){
                if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i_local ]=  (double) (-1);
                else phi[i_local ]= (PRECISION)(1);
                
                // debug print
                /*printf("local_i = %d on rank %d (i_global = %d)\n",
                      i_local , g.my_rank, i_global);*/
              }
            }
          }
        }
      }
      
    }else {
      error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
    }
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
  //exit(0);
}

/*for (int i_global = 0; i_global < N; i_global++) {
  owner = i_global / chunk;
  int i_local = i_global % chunk;
  if(g.my_rank==0) printf("component %d, by processor %d (local index %d)\n", i_global, owner, i_local);
}
  */  
/*void vector_PRECISION_define_random_rademacher( vector_PRECISION phi, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i]=  (double) (-1);
      else phi[i]= (PRECISION)(1);
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );
  
}
*/
