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


complex_PRECISION g5_3D_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = g5_3D_hutchinson_plain_PRECISION;
  estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
  trace += estimate.acc_trace/estimate.sample_size;

  // if deflation vectors are available
  //if(g.trace_deflation_type[l->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( l, h, threading );
  //}

  return trace;
}


complex_PRECISION g5_3D_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    } else {
      //vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    /* Apply Gamma_5 */
    gamma5_PRECISION( p->b, p->b, l, threading );

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
     // if(g.trace_deflation_type[l->depth] != 0){
     //   hutchinson_deflate_vector_PRECISION(p->x, l, threading);
     // }
      aux = global_inner_product_PRECISION( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );
    } else {
      //vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    }

    return aux;
  }
}

//-------------------------------- Jose’s code below -------------------------------------------------//

// if type_appl==-1 : Hutchinson-like
// else : direct term, where type_appl is the index of the deflation vector to apply
//        the operator on
complex_PRECISION hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
  
   // if ( type_appl==-1 ) {
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

    return aux;  
  }
}




// TODO: Unify with original mlmc driver (maybe by passing finest level l to blind)
complex_PRECISION mlmc_hutchinson_g5_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1; i++ ){ 
    // set the pointer to the mlmc difference operator
    h->hutch_compute_one_sample = hutchinson_mlmc_difference_PRECISION;
    estimate = hutchinson_blind_g5_PRECISION( l, lx->depth, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;
    lx = lx->next_level;
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  estimate = hutchinson_blind_g5_PRECISION( l, lx->depth, h, 0, threading );
  trace += estimate.acc_trace/estimate.sample_size;

  return trace;
}

// the term tr( A_{l}^{-1} - P A_{l+1}^{-1} R )
complex_PRECISION hutchinson_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // FIRST TERM : result stored in p->x
  // apply A_{l}^{-1}
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
  //if(l->depth == 1){printf("%d \n",l->depth);fflush(0);exit(0);}
    //TODO: Set a 3d-trace parameter in .ini --------------
    int d3 = 1;
    if(d3 == 1){
      vector_PRECISION_copy(h->mlmc_b1, h->rademacher_vector, start, end, l );  
      if( l->depth==0 ){
        gamma5_PRECISION( h->rademacher_vector, h->rademacher_vector, l, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, l->inner_vector_size, &startg5, &endg5, l, threading, l->num_lattice_site_var );
        coarse_gamma5_PRECISION( h->rademacher_vector, h->rademacher_vector, startg5, endg5, l );
      }
    }
    
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    } 
    // solution of this solve is in l->p_PRECISION.x
    apply_solver_PRECISION( l, threading );
  }

  // SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict
  // 2. invert
  // 3. Prolongate
  {
    //TODO: Set a 3d-trace parameter in .ini --------------    
    int d3 = 1;
    int start, end;
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    if(d3 == 1){//restore copy
      vector_PRECISION_copy(h->rademacher_vector, h->mlmc_b1, start, end, l );  
    }
    if ( type_appl==-1 ) {
      apply_R_PRECISION( l->next_level->p_PRECISION.b, h->rademacher_vector, l, threading );
    } 
    
    //TODO: Set a 3d-trace parameter in .ini --------------
    if(d3 == 1){  
      if( l->next_level->depth==0 ){
        gamma5_PRECISION( l->next_level->p_PRECISION.b, l->next_level->p_PRECISION.b, l->next_level, threading );
        //gamma5_PRECISION( l->next_level->p_PRECISION.b, l->next_level->p_PRECISION.b, l->next_level, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, l->next_level->inner_vector_size, &startg5, &endg5, l->next_level, threading, l->next_level->num_lattice_site_var );
        coarse_gamma5_PRECISION( l->next_level->p_PRECISION.b, l->next_level->p_PRECISION.b, startg5, endg5, l->next_level );
        //coarse_gamma5_PRECISION( l->next_level->p_PRECISION.b, l->next_level->p_PRECISION.b, startg5, endg5, l->next_level );
      }
    }
//exit(0);

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );
    apply_P_PRECISION( h->mlmc_b2, l->next_level->p_PRECISION.x, l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l);
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    
    int d3 = 1;
    if(d3 == 1){
      vector_PRECISION_copy(h->rademacher_vector, h->mlmc_b1, start, end, l );  
    }
    
    vector_PRECISION_minus( h->mlmc_b1, p->x, h->mlmc_b2, start, end, l); 

    if ( type_appl==-1 ) {
      aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l, threading );
    } 

    return aux; 
  }
}


struct sample hutchinson_blind_g5_PRECISION( level_struct *l, int depth, hutchinson_PRECISION_struct* h, int type, struct Thread *threading ){
  int i, j;
  complex_PRECISION one_sample=0, variance=0, trace=0;
  double RMSD;
  struct sample estimate;

  estimate.acc_trace = 0.0;

  int start, end;
  level_struct* lx = l;
  //get_correct_l_PRECISION( depth, lx );
  for (int d = 0; d < depth; d++){
    lx = lx->next_level;
  }
  //if(g.my_rank==0)printf("lx_depth = %d\n", lx->depth);
  // TODO : move this allocation to some init function
  complex_PRECISION* samples = (complex_PRECISION*) malloc( h->max_iters[lx->depth]*sizeof(complex_PRECISION) );
  memset( samples, 0.0, h->max_iters[lx->depth]*sizeof(complex_PRECISION) );

  level_struct* l_restrict = l;

  double t0 = MPI_Wtime();
  for( i=0; i<h->max_iters[lx->depth];i++ ){
    
// 1. create Rademacher vector, stored in h->rademacher_vector
    rademacher_create_PRECISION( l, h, type, threading );
    //compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    //vector_PRECISION_copy( h->mlmc_testing, h->rademacher_vector, start, end, l );
    l_restrict = l;  
    for (int d = 0; d < depth; d++){
       apply_R_PRECISION( h->rademacher_vector, h->rademacher_vector, l_restrict, threading );
       l_restrict = l_restrict->next_level;
       //if(g.my_rank==0)printf("restricting = %d times\n", d+1);fflush(0);
    /*  if(depth==1){
  printf("%d d=%d\n",l_restrict->depth,d);fflush(0);
     //exit(0);
  }*/
    //compute_core_start_end( 0, lx->inner_vector_size, &start, &end, lx, threading );
    //vector_PRECISION_copy(  h->rademacher_vector, h->mlmc_testing, start, end, lx );


    }

    // 2. apply the operator to the Rademacher vector
    // 3. dot product
    one_sample = h->hutch_compute_one_sample( -1, lx, h, threading );


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
      if( i > h->min_iters[lx->depth] && RMSD < cabs(trace) * h->trace_tol * h->tol_per_level[lx->depth]) break; 
    }
  }
  double t1 = MPI_Wtime();
  if(g.my_rank==0){
    printf("\n");
    printf("Time for sample computation (Avg.): \t %f\n\n", (t1-t0)/h->max_iters[lx->depth]);
  }
  estimate.sample_size = i;

  free(samples);

  
  return estimate;
}



// apply the interpolation
void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){
  if( l->depth==0 ){
    interpolate3_PRECISION( l->sbuf_PRECISION[0], in, l, threading );
    trans_back_PRECISION( (vector_double)out, l->sbuf_PRECISION[0], l->s_PRECISION.op.translation_table, l, threading );
  }
  else{
    interpolate3_PRECISION( out, in, l, threading );
  }
}


// apply the restriction
void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){
  if( l->depth==0 ){
    trans_PRECISION( l->sbuf_PRECISION[0], (vector_double)in, l->s_PRECISION.op.translation_table, l, threading );     
    restrict_PRECISION( out, l->sbuf_PRECISION[0], l, threading );
  }
  else{
    restrict_PRECISION( out, in, l, threading );
  }
}
