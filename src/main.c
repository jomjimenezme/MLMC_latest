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

#include "profiling.h"
#include "proxies/dirac_proxy.h"

#ifdef HAVE_HDF5
Hdf5_fileinfo h5info;
#endif
struct common_thread_data *commonthreaddata;

int main( int argc, char **argv ) {

  RangeHandleType profilingRangeMain = startProfilingRange("main");
    
#ifdef HAVE_HDF5
  h5info.filename=NULL;
  h5info.file_id=-1; 
  h5info.rootgroup_id=-1; 
  h5info.configgroup_id=-1;
  h5info.eigenmodegroup_id=-1;
  h5info.thiseigenmodegroup_id=-1;
  h5info.isOpen=0;
  h5info.mode=-1;
#endif
  level_struct l;
  config_double hopp = NULL, clov = NULL;
  
  MPI_Init( &argc, &argv );
  
  predefine_rank();
  if ( g.my_rank == 0 ) {
    printf("\n\n+------------------------------------------------------------------------+\n");
    printf("| The DDalphaAMG solver library.                                         |\n");
    printf("| Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Gustavo Ramirez, |\n");
    printf("|       Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori, Tilmann Matthaei, Ke-Long Zhang.   |\n");
    printf("|                                                                        |\n");
    printf("| This program comes with ABSOLUTELY NO WARRANTY.                        |\n");
    printf("+------------------------------------------------------------------------+\n\n");
  }

  method_init( &argc, &argv, &l );

  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);
  
  MALLOC( hopp, complex_double, 3*l.inner_vector_size );
  if ( g.two_cnfgs ) {
    MALLOC( clov, complex_double, 3*l.inner_vector_size );
    printf0("clover term configuration: %s", g.in_clov ); 

    if(g.in_format == _LIME)
      lime_read_conf( (double*)(clov), g.in_clov, &(g.plaq_clov) );
    else if(g.in_format == _MULTI)
      read_conf_multi( (double*)(clov), g.in, &(g.plaq_hopp), &l );
    else
      read_conf( (double*)(clov), g.in_clov, &(g.plaq_clov), &l );

    printf0("hopping term ");
  }

  if(g.in_format == _LIME)
    lime_read_conf( (double*)(hopp), g.in, &(g.plaq_hopp) );
  else if(g.in_format == _MULTI)
    read_conf_multi( (double*)(hopp), g.in, &(g.plaq_hopp), &l );
  else
    read_conf( (double*)(hopp), g.in, &(g.plaq_hopp), &l );

  if ( !g.two_cnfgs ) {
    g.plaq_clov = g.plaq_clov;
  }

  // store configuration, compute clover term
  dirac_setup( hopp, clov, &l );
  FREE( hopp, complex_double, 3*l.inner_vector_size );
  if ( g.two_cnfgs ) {
    FREE( clov, complex_double, 3*l.inner_vector_size );
  }

  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);

#pragma omp parallel num_threads(g.num_openmp_processes)
  {
    g.on_solve=0;

    struct Thread threading;
    l.threading = &threading;
    setup_threading(&threading, commonthreaddata, &l);
    setup_no_threading(no_threading, &l);

    // TODO: move this line to a better place !
    g.nr_threads = threading.n_core;
    RangeHandleType rangeHandle;
    
    rangeHandle = startProfilingRange("Setup");
    // setup up initial MG hierarchy
    method_setup( NULL, &l, &threading );
    endProfilingRange(rangeHandle);

    set_some_coarsest_level_improvs_params_for_setup( &l, &threading );

    rangeHandle = startProfilingRange("Update");
    // iterative phase
    method_update( l.setup_iter, &l, &threading );
    endProfilingRange(rangeHandle);

    set_some_coarsest_level_improvs_params_for_solve( &l, &threading );

    rangeHandle = startProfilingRange("Solve");
    g.on_solve=1;

    //solve_driver( &l, &threading );
    // pre-setting some values for the Hutchinson struct
    struct Thread *threadingx = &threading;  
    complex_double trace;
    int rnd_seed = 1234;
    srand( time( 0 ) + rnd_seed*g.my_rank );
    {  
      hutchinson_diver_double_init( &l, &threading );  
      hutchinson_diver_double_alloc( &l, &threading ); 
    }
    for(g.time_slice = 0; g.time_slice <4 ; g.time_slice++){//g.time_slice< g.global_lattice[0][0]; g.time_slice++){
      
      if(g.my_rank==0) printf("\n\n Timeslice %d\n\n",  g.time_slice);
      
     if( g.trace_op_type == 2){
        START_MASTER(threadingx)
        if(g.my_rank==0) printf("Using plain Hutchinson for computing the trace\n");
        END_MASTER(threadingx)
          
        trace = hutchinson_driver_double( &l, &threading );
        
        START_MASTER(threadingx)
        if(g.my_rank==0) printf("\nResulting trace from plain Hutchinson = %f+i%f\n", CSPLIT(trace)); fflush(0); 
        END_MASTER(threadingx)
      }
      
      if( g.trace_op_type == 1){

        START_MASTER(threadingx)
        if(g.my_rank==0) printf("Using MLMC for computing the trace\n");
        END_MASTER(threadingx)
          
        trace = mlmc_hutchinson_g5_driver_double(&l, &threading);
          
        START_MASTER(threadingx)
        if(g.my_rank==0) printf("\nResulting trace from MLMC  = %f+i%f\n", CSPLIT(trace)); fflush(0); 
        END_MASTER(threadingx)
      }

    }
    
    hutchinson_diver_double_free( &l, &threading );  
    
    endProfilingRange(rangeHandle);
  }

  finalize_common_thread_data(commonthreaddata);
  finalize_no_threading(no_threading);

  method_free( &l );
  method_finalize( &l );

  MPI_Finalize();
  
  endProfilingRange(profilingRangeMain);
  return 0;
}
