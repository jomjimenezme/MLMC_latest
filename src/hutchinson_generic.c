#include "main.h"

struct estimate {
  int counter; //required number of estimates.
  complex_PRECISION estimate;
};

struct sample {
  // required number of estimates
  int sample_size;
  // accumulated trace
  complex_PRECISION acc_trace;
};

void hutchinson_diver_PRECISION_init( level_struct *l, struct Thread *threading ) {
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);

  h->max_iters = NULL;
  h->min_iters = NULL;
        
  // MLMC
  h->mlmc_b1 = NULL;
  h->mlmc_b2 =NULL;
  h->mlmc_testing =NULL;
  h->rademacher_vector =NULL;
  SYNC_MASTER_TO_ALL(threading)
}

void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading ) {
  int i;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;

  PUBLIC_MALLOC( h->max_iters, int, g.num_levels );
  PUBLIC_MALLOC( h->min_iters, int, g.num_levels );

  //For MLMC
  PUBLIC_MALLOC( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->mlmc_testing, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );

  for ( i=0;i<g.num_levels;i++ ) {
    h->max_iters[i] = g.trace_max_iters[i];
    h->min_iters[i] = g.trace_min_iters[i];
  }
}

void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading ) {
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;

  PUBLIC_FREE( h->max_iters, int, g.num_levels );
  PUBLIC_FREE( h->min_iters, int, g.num_levels );
        
  PUBLIC_FREE( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );   
  PUBLIC_FREE( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );   
  PUBLIC_FREE( h->mlmc_testing, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );   
}

complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){    
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
  trace += estimate.acc_trace/estimate.sample_size;

  
  return trace;
}

void rademacher_create_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading ){
  if( type==0 ){
    START_MASTER(threading)
    vector_PRECISION_define_random_rademacher( h->rademacher_vector, 0, l->inner_vector_size, l );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
  }
  else if( type==1 ){
    START_MASTER(threading)
    vector_PRECISION_define_random_rademacher( h->rademacher_vector, 0, l->next_level->inner_vector_size, l->next_level );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
   }
  else{ error0("Unknown value for type of Rademacher vector in relation to level of creation\n"); }
}

// FIXME : hacked this function a bit to avoid compiler warnings
gmres_PRECISION_struct* get_p_struct_PRECISION( level_struct* l ){
  if( l->depth==0 ){
    if ( strcmp("PRECISION","float")==0 ) { return &(l->p_PRECISION); }
    else { return (gmres_PRECISION_struct*)&(g.p); }
  }
  else{ return &(l->p_PRECISION); }
}

int apply_solver_PRECISION( level_struct* l, struct Thread *threading ){
  int nr_iters;
  double buff1=0, buff2=0;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    
  buff1 = p->tol;
  p->tol = g.tol;
  START_MASTER(threading);
  if( l->level==0 ){
    buff2 = g.coarse_tol;
    g.coarse_tol = g.tol;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);
    
  nr_iters = fgmres_PRECISION( p, l, threading );

  START_MASTER(threading);
  p->tol = buff1;
  if( l->level==0 ){
    g.coarse_tol = buff2;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);

  return nr_iters;
}

// type : in case of 0 create Rademacher vectors at level l, in case of 1 create Rademacher vectors at level l->next_level
struct sample hutchinson_blind_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading ){
  int i, j;
  complex_PRECISION one_sample=0.0, variance=0.0, trace=0.0;
  double RMSD;
  struct sample estimate;

  // TODO : move this allocation to some init function
  complex_PRECISION* samples = (complex_PRECISION*) malloc( h->max_iters[l->depth]*sizeof(complex_PRECISION) );
  memset( samples, 0.0, h->max_iters[l->depth]*sizeof(complex_PRECISION) );

  estimate.acc_trace = 0.0;
  double t0 = MPI_Wtime();
  
  for( i=0; i<h->max_iters[l->depth];i++ ){
    // 1. create Rademacher vector, stored in h->rademacher_vector
    rademacher_create_PRECISION( l, h, type, threading );

    // 2. apply the operator to the Rademacher vector
    // 3. dot product
    one_sample = h->hutch_compute_one_sample( -1, l, h, threading );

    samples[i] = one_sample;

    // 4. compute estimated trace and variance, print something?
    estimate.acc_trace += one_sample;

    if( i!=0 ){
      variance = 0.0;
      estimate.sample_size = i+1;
      trace = estimate.acc_trace/estimate.sample_size;
      for( j=0; j<i; j++ ){
        variance += conj(samples[j] - trace) * (samples[j] - trace);
      }
      variance = variance / j;
      START_MASTER(threading);
      if(g.my_rank==0) {
        printf("[%d, trace: %f+%f, variance: %f] ", i, creal(trace), cimag(trace), creal(variance));
        fflush(0);
      }
      END_MASTER(threading);
      RMSD = sqrt(creal(variance)/j);
      if( i > h->min_iters[l->depth] && RMSD < cabs(trace) * h->trace_tol * h->tol_per_level[l->depth]) break; 
    }
  }
  double t1 = MPI_Wtime();
  if(g.my_rank==0){
    printf("\n");
    printf("Time for sample computation (Avg.): \t %f\n\n", (t1-t0)/h->max_iters[l->depth]);
  }

  estimate.sample_size = i;

  free(samples);

  return estimate;
}




// if type_appl==-1 : Hutchinson-like
// else : direct term, where type_appl is the index of the deflation vector to apply
//        the operator on
complex_PRECISION hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    PRECISION norm = 0.0;
   /* vector_PRECISION_define( p->b, 1.0, start, end, l );
    norm = global_norm_PRECISION( p->b, 0, l->inner_vector_size, l, threading );
    if(g.my_rank==0)
      printf("\n\t  Test 1: ||b|| = %e\n", norm);

    apply_solver_PRECISION( l, threading );

    norm = global_norm_PRECISION( p->x, 0, l->inner_vector_size, l, threading );

    if(g.my_rank==0)
      printf("\n\t  Test 1: ||x|| = %e\n", norm);
    */
    //------------------------------------------------------
    vector_PRECISION_define_spin_color( p->b, 0, l->inner_vector_size, l , threading);
    norm = global_norm_PRECISION( p->b, 0, l->inner_vector_size, l, threading );
    if(g.my_rank==0){
	 START_LOCKED_MASTER(threading)
      printf("\n\t  Test 2: ||b|| = %e\n", norm);
    END_LOCKED_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
   }

    apply_solver_PRECISION( l, threading );
    norm = global_norm_PRECISION( p->x, 0, l->inner_vector_size, l, threading );
    if(g.my_rank==0){
      START_LOCKED_MASTER(threading)

      printf("\n\t  Test 2: ||x|| = %e\n", norm);
      END_LOCKED_MASTER(threading)
     SYNC_MASTER_TO_ALL(threading)
    }
    
    
    //------------------------------------------------------
    vector_PRECISION_define_spin_color( h->rademacher_vector, start, end, l, threading );
    apply_operator_PRECISION(  p->b, h->rademacher_vector, p, l, threading );
    //vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    
    norm = global_norm_PRECISION( p->b, 0, l->inner_vector_size, l, threading );
    if(g.my_rank==0){
      START_LOCKED_MASTER(threading)

      printf("\n\t  Test 3: ||b|| = %e\n", norm);
      END_LOCKED_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    
    apply_solver_PRECISION( l, threading );
    norm = global_norm_PRECISION( p->x, 0, l->inner_vector_size, l, threading );
    if(g.my_rank==0){
      START_LOCKED_MASTER(threading)

      printf("\n\t  Test 3: ||x|| = %e\n", norm);
        END_LOCKED_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }

    vector_PRECISION_minus(p->b , p->x, h->rademacher_vector, start, end, l );
    PRECISION error = global_norm_PRECISION( p->b, 0, l->inner_vector_size, l, threading );
    if(g.my_rank==0){
      START_LOCKED_MASTER(threading)
      printf("\n\t  Test 3:  ||error||/||x|| = %e\n", error/norm);
      END_LOCKED_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    
    

    

  }

  complex_PRECISION aux = 0.0;
   /*// if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    //} else {
     // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    //}
   
    int d3 = 1;
    if(d3 == 1){
      if( l->depth==0 ){
        gamma5_PRECISION( p->b, p->b, l, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, l->inner_vector_size, &startg5, &endg5, l, threading, l->num_lattice_site_var );
        coarse_gamma5_PRECISION( p->b, p->b, startg5, endg5, l );
	//coarse_gamma5_PRECISION( p->b, p->b, startg5, endg5, l );
      }
    } 
  }

  {
    apply_solver_PRECISION( l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    
   
    if ( type_appl==-1 ) {
//      if(g.trace_deflation_type[l->depth] != 0){
//        hutchinson_deflate_vector_PRECISION(p->x, l, threading); 
//      }
      aux = global_inner_product_PRECISION( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );  
      } 
//    } else {
//      vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);  
//      aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
//    }
  */
    return aux;  

}

