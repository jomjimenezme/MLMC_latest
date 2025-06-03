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
/*  Rademacher ±1 point-source on a fixed global time-slice
 *  -------------------------------------------------------
 *  • Layout   :  site-major,  dof_idx = spin*3 + colour
 *  • Traversal:  only the sites owned by this rank
 *                lexicographic (t,z,y,x) with x fastest
 *  • Value    :  ±1   (real)  with probability ½ each,   imag = 0
 *  • The source is non-zero only on the global time-slice   g.time_slice.
 */
void vector_PRECISION_define_random_rademacher(vector_PRECISION phi,
                                               int start, int end,
                                               level_struct *l)
{
	int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
   if(thread == 0){
    if ( phi != NULL ) {
      int depth = l->depth;
      int r_t = g.my_coords[T];
      int r_z = g.my_coords[Z];
      int r_y = g.my_coords[Y];
      int r_x = g.my_coords[X];
      //printf("rank= %d \t %d %d %d %d \n",g.my_rank, r_t, r_z, r_y, r_x);

      /* global_splitting: Num Procs in each dim
       * Nt_loc: Local Lattice sites per direcction*/
      int Nt_loc = g.global_lattice[depth][T]/l->global_splitting[T];
      int Nz_loc = g.global_lattice[depth][Z]/l->global_splitting[Z];
      int Ny_loc = g.global_lattice[depth][Y]/l->global_splitting[Y];
      int Nx_loc = g.global_lattice[depth][X]/l->global_splitting[X];

      // Initial global coordinates based on cartesian coordinates
      int t0 = r_t * Nt_loc;
      int x0 = r_x * Nx_loc;
      int y0 = r_y * Ny_loc;
      int z0 = r_z * Nz_loc;

      //printf("RANK = %d \t Nt= %d,  Nz = %d, Ny= %d, Nx = %d\n T = %d, Z = %d, Y = %d, X = %d\n", g.my_rank, Nt_loc, Nz_loc, Ny_loc, Nx_loc, T,Z,Y,X);
      for (int lt=0; lt<Nt_loc; ++lt) {
        int t = t0 + lt;               // global t
        bool ontimeslice = (t ==g.time_slice);
        for (int lz=0; lz<Nz_loc; ++lz) {
          int z = z0 + lz;
            for (int ly=0; ly<Ny_loc; ++ly) {
              int y = y0 + ly;           // global y
                for (int lx=0; lx<Nx_loc; ++lx) {
                  int x = x0 + lx;             // global x

                  int loc_site  = ((lt * Nz_loc + lz) * Ny_loc + ly) * Nx_loc + lx;

                  for (int d=0; d<4; ++d) {
                    for (int c=0; c<3; ++c){
                        PRECISION re =  t
                                      + 3.14 * z
                                      - 2.72 * y
                                      + 0.58 * x
                                      - 1.41 * c
                                      + 1.20 * d;

                        int i = loc_site*12    /* site base  */
                                  + d * 3      /* stride is 3 colors  */
                                  + c;         /* color component   */
                        PRECISION non_re = 1.0 / (1.2345 + re);
                        if (ontimeslice){
                        	if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i ]=  (double) (-1);
               						else phi[i ]= (PRECISION)(1);}
                        else{phi[i]=0.0;}
                }
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

}

void vector_PRECISION_define_spin_color( vector_PRECISION phi, int start, int end, level_struct *l, struct Thread *threading) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
   if(thread == 0){
    if ( phi != NULL ) {
      int depth = l->depth;
      int r_t = g.my_coords[T];
      int r_z = g.my_coords[Z];
      int r_y = g.my_coords[Y];
      int r_x = g.my_coords[X];
      //printf("rank= %d \t %d %d %d %d \n",g.my_rank, r_t, r_z, r_y, r_x);
      
      /* global_splitting: Num Procs in each dim
       * Nt_loc: Local Lattice sites per direcction*/      
      int Nt_loc = g.global_lattice[depth][T]/l->global_splitting[T]; 
      int Nz_loc = g.global_lattice[depth][Z]/l->global_splitting[Z];
      int Ny_loc = g.global_lattice[depth][Y]/l->global_splitting[Y];
      int Nx_loc = g.global_lattice[depth][X]/l->global_splitting[X];
      
      // Initial global coordinates based on cartesian coordinates
      int t0 = r_t * Nt_loc;
      int x0 = r_x * Nx_loc;
      int y0 = r_y * Ny_loc;
      int z0 = r_z * Nz_loc;    
      
      //printf("RANK = %d \t Nt= %d,  Nz = %d, Ny= %d, Nx = %d\n T = %d, Z = %d, Y = %d, X = %d\n", g.my_rank, Nt_loc, Nz_loc, Ny_loc, Nx_loc, T,Z,Y,X);
      for (int lt=0; lt<Nt_loc; ++lt) {
        int t = t0 + lt;               // global t
        for (int lz=0; lz<Nz_loc; ++lz) {
          int z = z0 + lz;
            for (int ly=0; ly<Ny_loc; ++ly) {
              int y = y0 + ly;           // global y
                for (int lx=0; lx<Nx_loc; ++lx) {
                  int x = x0 + lx;             // global x

                  int loc_site  = ((lt * Nz_loc + lz) * Ny_loc + ly) * Nx_loc + lx;

                  for (int d=0; d<4; ++d) {
                    for (int c=0; c<3; ++c){
                        PRECISION re =  t
                                      + 3.14 * z
                                      - 2.72 * y
                                      + 0.58 * x
                                      - 1.41 * c
                                      + 1.20 * d;

                        int i = loc_site*12    /* site base  */
                                  + d * 3      /* stride is 3 colors  */
                                  + c;         /* color component   */
                        PRECISION non_re = 1.0 / (1.2345 + re); 
                        phi[i] = re + I* non_re;;
                }
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



/*void vector_PRECISION_define_spin_color( vector_PRECISION phi, int start, int end, level_struct *l, struct Thread *threading) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
   if(thread == 0){
    if ( phi != NULL ) {
      int i_global, i_local, site_index, owner;
      int dof_per_site = l->num_lattice_site_var; // variables per site
      
      int depth = l->depth;
      int global_sites = g.global_lattice[depth][0] *
      g.global_lattice[depth][1] *
      g.global_lattice[depth][2] *
      g.global_lattice[depth][3];
      
      //  degrees of freedom in the full global vector
      int total_dof = dof_per_site * global_sites;
      
      // Number of dofs owned by each process //TODO: is always even?
      int chunk = total_dof / g.num_processes;
      
      for (int t = 0; t < g.global_lattice[0][0]; t++) {
        for (int x = 0; x < g.global_lattice[0][3]; x++) {
          for (int y = 0; y < g.global_lattice[0][2]; y++) {
            for (int z = 0; z < g.global_lattice[0][1]; z++) {
              
              // Get the global site index at (t_slice, z, y, x)
              site_index = lex_index(t, x, y, z, g.global_lattice[0]);
              if(g.my_rank==2)printf("[%d, %d, %d, %d]\t site_index \t\t %d\n",
                          t,z,y,x, site_index);
              // Loop over the dof of that site
              int counter = 0; int color = 0; int spin = 0;
              for (i_global = dof_per_site * site_index;
                  i_global < dof_per_site * site_index + dof_per_site;
                  i_global++) {
                  
                  if(color == 3){ 
                    color = 0;
                    spin ++;
                  }
                  // Compute owner of the entry and local index 
                  owner = i_global / chunk;
                  i_local = i_global % chunk;
                  
                  // Only the owner assigns the value
                  if(owner == g.my_rank){
//START_LOCKED_MASTER(threading);
                     PRECISION re = t + x*3.14 - y*2.72 + z*0.58 - color*1.41 + spin*1.20;
                     PRECISION non_re = 1.0 / (1.2345 + re);
                     phi[i_local ] = re + I* non_re;
                    // debug print
                    //printf("[%d, %d, %d, %d]\t c = %d s= %d \t\t %d\n",
                    //     t,z,y,x,color,spin,counter);
//END_LOCKED_MASTER(threading);
//SYNC_MASTER_TO_ALL(threading);
                  }
                  
                  color++;
                  counter++;
                
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
*/

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
