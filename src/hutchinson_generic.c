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

  h->finest_level = l;
  h->lx_i = NULL;
  h->lx_j = NULL;

  SYNC_MASTER_TO_ALL(threading)
}


void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading ) {
  int i;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;

  PUBLIC_MALLOC( h->max_iters, int, g.num_levels );
  PUBLIC_MALLOC( h->min_iters, int, g.num_levels );

  // For MLMC
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


/*

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
*/


void rademacher_create_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading ){
  if( type==0 ){
    START_MASTER(threading)
    if(h->hutch_compute_one_sample == g5_3D_hutchinson_mlmc_difference_PRECISION || h->hutch_compute_one_sample == g5_3D_hutchinson_mlmc_coarsest_PRECISION){
      vector_PRECISION_define_random_rademacher( h->rademacher_vector, 0, h->finest_level->inner_vector_size, h->finest_level );
    }else{
      vector_PRECISION_define_random_rademacher( h->rademacher_vector, 0, l->inner_vector_size, l );
    }
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
  int nr_iters = 0;
  double buff1=0, buff2=0;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

  p->print_iters = 1;

  buff1 = p->tol;
  p->tol = g.tol;
  START_MASTER(threading);
  if( l->level==0 ){
    buff2 = g.coarse_tol;
    g.coarse_tol = g.tol;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);

  if ( l->level > 0 ) {
    nr_iters = fgmres_PRECISION( p, l, threading );
  } else {
#ifdef GCRODR
    // NOTE : something that shouldn't be happening here happens, namely the RHS is changed
    //        by the function coarse_solve_odd_even_PRECISION(...). So, we back it up and restore
    //        it as necessary

    int start,end;
    //compute_core_start_end( l->p_PRECISION.v_start, l->p_PRECISION.v_end, &start, &end, l, threading );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( l->p_PRECISION.rhs_bk, l->p_PRECISION.b, start, end, l );

    START_MASTER(threading)
    l->p_PRECISION.was_there_stagnation = 0;
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    while( 1 ) {
      coarse_solve_odd_even_PRECISION( &(l->p_PRECISION), &(l->oe_op_PRECISION), l, threading );
      if ( l->p_PRECISION.was_there_stagnation==0 ) { break; }
      else if ( l->p_PRECISION.was_there_stagnation==1 && l->p_PRECISION.gcrodr_PRECISION.CU_usable==1 ) {
        // in case there was stagnation, we need to rebuild the coarsest-level data
        double time_bk = g.coarsest_time;
        coarsest_level_resets_PRECISION( l, threading );
        START_MASTER(threading)
        l->p_PRECISION.was_there_stagnation = 0;
        g.coarsest_time = time_bk;
        END_MASTER(threading)
        SYNC_MASTER_TO_ALL(threading)
        vector_PRECISION_copy( l->p_PRECISION.b, l->p_PRECISION.rhs_bk, start, end, l );
      }
      else {
        // in this case, there was stagnation but no deflation/recycling subspace is being used
        break;
      }
    }
#else
    coarse_solve_odd_even_PRECISION( &(l->p_PRECISION), &(l->oe_op_PRECISION), l, threading );
#endif
  }

  START_MASTER(threading);
  p->tol = buff1;
  if( l->level==0 ){
    g.coarse_tol = buff2;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);

  p->print_iters = 0;

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

      if(i == h->max_iters[l->depth] - 1 && g.trace_op_type != 7)
        g.variances[l->depth] += creal(variance);

      if(i == h->max_iters[l->depth] - 1 && g.trace_op_type == 7){
	int nlevs = g.num_levels;
	int idx = h->lx_i->depth*nlevs + h->lx_j->depth;
	g.variances[idx] += creal(variance);
	//printf("\nOperator index in variance assignment = %d", idx);
      }
      }
      END_MASTER(threading);
      RMSD = sqrt(creal(variance)/j);
      if( i > h->min_iters[l->depth] && RMSD < cabs(trace) * h->trace_tol * h->tol_per_level[l->depth]) break;
    }
  }
  double t1 = MPI_Wtime();
  if(g.my_rank==0) {
    printf("\n");
    printf("Time for sample computation (Avg.): \t %f\n\n", (t1-t0)/h->max_iters[l->depth]);
  }

  estimate.sample_size = i;

  free(samples);

  return estimate;
}

// this is the driver for plain Hutchinson
complex_PRECISION g5_3D_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = g5_3D_hutchinson_plain_PRECISION;

  if (g.probing) {
    for (g.coloring_count = 1; g.coloring_count < g.num_colors[0] + 1; g.coloring_count++){
      if(g.my_rank==0) {
        printf("\nEstimating trace at color %d\n", g.coloring_count);
      }

      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
      trace += estimate.acc_trace / estimate.sample_size;
    }
  } else {
    estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
    trace += estimate.acc_trace / estimate.sample_size;
  }

  // if deflation vectors are available
  //if(g.trace_deflation_type[l->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( l, h, threading );
  //}

  return trace;
}

// this is the driver for plain Hutchinson, BUT for the connected operator (-> connected diagram)
complex_PRECISION g5_3D_connected_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = g5_3D_connected_hutchinson_plain_PRECISION;
  for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected<g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {
    
    if (g.probing) {
      for (g.coloring_count = 1; g.coloring_count < g.num_colors[0] + 1; g.coloring_count++){
        if(g.my_rank==0) {
          printf("\nEstimating trace at color %d\n", g.coloring_count);
        }

        estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
        trace += estimate.acc_trace / estimate.sample_size;
      }
    } else {
      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
      trace += estimate.acc_trace / estimate.sample_size;
    }
    // estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
    // trace += estimate.acc_trace/estimate.sample_size;
  }

  // if deflation vectors are available
  //if(g.trace_deflation_type[l->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( l, h, threading );
  //}

  return trace;
}


complex_PRECISION g5_3D_connected_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      //Apply \Pi_{t’+t}
      int bufft = g.time_slice;
      // TODO : check : is this assuming periodic in time ?
      g.time_slice = g.time_slice + g.time_slice_inner_connected;
      g.time_slice = g.time_slice%g.global_lattice[0][0];
      vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
      g.time_slice = bufft;
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    } else {
      //vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Apply Gamma_5
    gamma5_PRECISION( p->b, p->b, l, threading );
  }

  {
    apply_solver_PRECISION( l, threading );
  }

  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    // apply \Pi_{t'}
    //vector_PRECISION_ghg(  p->x, 0, l->inner_vector_size, l );
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    //g.time_slice = g.time_slice%g.global_lattice[0][0];
    //vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
    vector_PRECISION_ghg( p->x, 0, l->inner_vector_size, l );
    g.time_slice = bufft;

    // apply G5
    gamma5_PRECISION( p->x, p->x, l, threading );

    // inverse again
    vector_PRECISION_copy( p->b, p->x, start, end, l );
    apply_solver_PRECISION( l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux = 0;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
    //if(g.trace_deflation_type[l->depth] != 0){
    //  hutchinson_deflate_vector_PRECISION(p->x, l, threading);
    //}
      aux = global_inner_product_PRECISION( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );
    } else {
      //vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    }

    return aux;
  }
}

complex_PRECISION g5_3D_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      //vector_PRECISION_copy( h->mlmc_b1, h->rademacher_vector, start, end, l );
      vector_PRECISION_ghg(  h->rademacher_vector, 0, l->inner_vector_size, l );
      vector_PRECISION_copy( p->b,  h->rademacher_vector, start, end, l );
      //vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    } else {
      //vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Apply Gamma_5
    gamma5_PRECISION( p->b, p->b, l, threading );
  }

  {
    apply_solver_PRECISION( l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux = 0;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
    //if(g.trace_deflation_type[l->depth] != 0){
    //  hutchinson_deflate_vector_PRECISION(p->x, l, threading);
    //}
      aux = global_inner_product_PRECISION( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );
    } else {
      //vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    }

    return aux;
  }
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





























complex_PRECISION g5_3D_hutchinson_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }

  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );

  // FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // Project rademacher vector into time-slice
    vector_PRECISION_ghg(  h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
      
    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level

    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }

    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }

  // SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict (iteratively from finest)
  // 2. invert
  // 3. Prolongate (iteratively to finest)
  {
    int start, end, start_tmp, end_tmp, start_tmp_c, end_tmp_c;
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    if ( type_appl==-1 ) {
      // Copy gamma_5 x at the finest
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

      // Restrict it from finest to NEXT level
      for (int i = 0; i < l_index+1; ++i) {
        level_struct *fine   = levels[i];
        level_struct *coarse = levels[i + 1];

        compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
        compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
        // restrict from fine to coarse
        apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
        vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

        if (g.my_rank == 0)
          printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
      }

      // Copy restricted vector to rhs at the next level
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
      vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, end_tmp, l->next_level);
    } else {
    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );

    // Prolongate Solution from NEXT level to finest

    vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
    for (int i = l_index+1; i > 0; --i) {
      level_struct *coarse = levels[i];
      // finer level (towards finest)
      level_struct *fine   = levels[i - 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // prolongate from coarse to fine
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      // copy result to fine level
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

      if(g.my_rank == 0) {
        printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in second term vector
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    //copy for dot product
    vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }

  // subtract the results and perform dot product AT FINEST
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, finest_l);

    if ( type_appl==-1 ) {
      //if(g.trace_deflation_type[l->depth] != 0){
        //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
      //}
      aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, finest_l, threading );
    } else {
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
  }

  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );

  return aux;
}


complex_PRECISION g5_3D_connected_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }

  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );

  // CONNECTED STEP #1. apply \Pi_{t+t'}
  // CONNECTED STEP #2. apply G5

  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
    //vector_PRECISION_ghg(  h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );
  }

  // CONNECTED STEP #3. apply level difference operator

  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level
    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // LEVEL DIFFERENCE : SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict (iteratively from finest)
  // 2. invert
  // 3. Prolongate (iteratively to finest)
  {
    int start, end, start_tmp, end_tmp, start_tmp_c, end_tmp_c;
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    if ( type_appl==-1 ) {
      // Copy gamma_5 x at the finest
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

      // Restrict it from finest to NEXT level
      for (int i = 0; i < l_index+1; ++i) {
        level_struct *fine   = levels[i];
        level_struct *coarse = levels[i + 1];

        compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
        compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
        // restrict from fine to coarse
        apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
        vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

        if (g.my_rank == 0)
          printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
      }

      // Copy restricted vector to rhs at the next level
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
      vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, end_tmp, l->next_level);
    } else {
    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );

    // Prolongate Solution from NEXT level to finest
    vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
    for (int i = l_index+1; i > 0; --i) {
      level_struct *coarse = levels[i];
      // finer level (towards finest)
      level_struct *fine   = levels[i - 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // prolongate from coarse to fine
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      // copy result to fine level
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

      if(g.my_rank == 0) {
        printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in second term vector
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    //copy for dot product
    vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // subtract the results
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, finest_l);
  }

  // CONNECTED STEP #4. apply \Pi_{t'}
  // CONNECTED STEP #5. apply G5

  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // apply \Pi_{t'}
    //vector_PRECISION_ghg(  p->x, 0, l->inner_vector_size, l );
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    //g.time_slice = g.time_slice%g.global_lattice[0][0];
    //vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
    vector_PRECISION_ghg( h->mlmc_b1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->mlmc_b1, finest_l, threading );
  }

  // CONNECTED STEP #6. apply level difference operator

  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level

    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }

    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // LEVEL DIFFERENCE : SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict (iteratively from finest)
  // 2. invert
  // 3. Prolongate (iteratively to finest)
  {
    int start, end, start_tmp, end_tmp, start_tmp_c, end_tmp_c;
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    if ( type_appl==-1 ) {
      // Copy gamma_5 x at the finest
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

      // Restrict it from finest to NEXT level
      for (int i = 0; i < l_index+1; ++i) {
        level_struct *fine   = levels[i];
        level_struct *coarse = levels[i + 1];

        compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
        compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
        // restrict from fine to coarse
        apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
        vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

        if (g.my_rank == 0)
          printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
      }

      // Copy restricted vector to rhs at the next level
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
      vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, end_tmp, l->next_level);
    } else {
    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );

    // Prolongate Solution from NEXT level to finest

    vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
    for (int i = l_index+1; i > 0; --i) {
      level_struct *coarse = levels[i];
      // finer level (towards finest)
      level_struct *fine   = levels[i - 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // prolongate from coarse to fine
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      // copy result to fine level
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

      if(g.my_rank == 0) {
        printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in second term vector
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    //copy for dot product
    vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // subtract the results
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, finest_l);
  }

  // do the dot product AT FINEST
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);

    if ( type_appl==-1 ) {
      //if(g.trace_deflation_type[l->depth] != 0){
        //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
      //}
      aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, finest_l, threading );
    } else {
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
  }

  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );

  return aux;
}


complex_PRECISION g5_3D_hutchinson_mlmc_coarsest_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }

  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );

  // FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // Project rademacher vector into time-slice
    vector_PRECISION_ghg(  h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level

    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // TODO : CHECK if coarse or fine !!
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }

  // // SECOND TERM : result stored in h->mlmc_b2
  // // 1. Restrict (iteratively from finest)
  // // 2. invert
  // // 3. Prolongate (iteratively to finest)
  // {
  //   int start, end, start_tmp, end_tmp;
  //   compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

  //   if ( type_appl==-1 ) {
  //     // Copy gamma_5 x at the finest
  //     vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

  //     // Restrict it from finest to NEXT level
  //     for (int i = 0; i < l_index+1; ++i) {
  //       level_struct *fine   = levels[i];
  //       level_struct *coarse = levels[i + 1];

  //       compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
  //       compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
  //       // restrict from fine to coarse
  //       apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
  //       vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

  //       if (g.my_rank == 0)
  //         printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
  //     }

  //     // Copy restricted vector to rhs at the next level
  //     compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
  //     vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, start_tmp, l->next_level);
  //   } else {
  //    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
  //   }

  //   // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
  //   apply_solver_PRECISION( l->next_level, threading );

  //   // Prolongate Solution from NEXT level to finest

  //   vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
  //   for (int i = l_index+1; i > 0; --i) {
  //     level_struct *coarse = levels[i];
  //     // finer level (towards finest)
  //     level_struct *fine   = levels[i - 1];

  //     compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
  //     // prolongate from coarse to fine
  //     apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
  //     // copy result to fine level
  //     vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

  //     if(g.my_rank == 0)
  //         printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
  //   }

  //   // Copy result in second term vector
  //   compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
  //   //copy for dot product
  //   vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  // }

  // subtract the results and perform dot product AT FINEST
  {
    // int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    // compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    // vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, finest_l);

    if ( type_appl==-1 ) {
      //if(g.trace_deflation_type[l->depth] != 0){
        //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
      //}
      aux = global_inner_product_PRECISION( h->rademacher_vector, first_term, p->v_start, p->v_end, finest_l, threading );
    } else {
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
  }

  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );

  return aux;
}


// complex_PRECISION g5_3D_hutchinson_mlmc_coarsest_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
//   /* Store from fine to coarsest lvl in an array*/
//   level_struct *finest_l = h->finest_level;
//   level_struct *levels[g.num_levels];
//   int lvl_nr = 0, l_index = -1;

//   if(g.my_rank == 0)
//     printf("\n\n------------------function at COARSEST%d ------------------\n\n", l->depth);fflush(0);

//   for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
//     levels[lvl_nr] = l_tmp;
//     if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
//     if (l_tmp == l) l_index = lvl_nr;
//     lvl_nr++;
//   }

//  // Preallocate coarse term vector at the finest
//   vector_PRECISION *only_term;
//   PUBLIC_MALLOC( only_term, complex_PRECISION, finest_l->inner_vector_size );

//   // FIRST TERM : result stored in only_term
//   {
//     int start, end;
//     //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l );
//     compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

//     /* Apply gamma_5 AT FINEST level and save in buffer */
//     gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );

//     /* copy buffer to be restricted below */
//     if ( type_appl==-1 ) {
//       vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
//     } else {
//      // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
//     }

//     /* Restrict rhs from finest to current level */
//     int start_tmp, end_tmp;

//     for (int i = 0; i < l_index; ++i) {
//       level_struct *fine   = levels[i];
//       level_struct *coarse = levels[i + 1];

//       compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
//       apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  // restrict from fine to coarse
//       vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, coarse);

//       if (g.my_rank == 0)
//         printf("Coarsest 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
//     }
//     /* Solve at current level */
//     compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
//     gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
//     vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);

//     apply_solver_PRECISION( l, threading );

//     /* Prolongate solution to finest level */
//     vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);

//     for (int i = l_index; i > 0; --i) {
//       level_struct *coarse = levels[i];
//       level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

//       compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
//       apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
//       vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

//       if(g.my_rank == 0)
//           printf("Coarsest: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);

//     }

//     // Copy result in the term vector
//     compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
//     vector_PRECISION_copy(only_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);  /* copy for dot product */

//   }

//   //  perform dot product with rademacher AT FINEST
//   {
//     int start, end;
//     complex_PRECISION aux;
//     gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l );

//     if ( type_appl==-1 ) {
//       //if(g.trace_deflation_type[l->depth] != 0){
//         //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
//       //}
//       aux = global_inner_product_PRECISION( h->rademacher_vector, only_term, p->v_start, p->v_end, finest_l, threading );
//     } else {
//       //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
//     }

//     return aux;
//   }

// }


complex_PRECISION g5_3D_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1; i++ ){
    // set the pointer to the mlmc difference operator
    h->hutch_compute_one_sample = g5_3D_hutchinson_mlmc_difference_PRECISION;
    if (g.probing) {
      for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++){
        if(g.my_rank==0)
          printf("\nEstimating trace at color %d\n", g.coloring_count);
      
        estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
        trace += estimate.acc_trace / estimate.sample_size;
      }
    } else {
      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
      trace += estimate.acc_trace / estimate.sample_size;
    }
    
    //if deflation vectors are available
    //if(g.trace_deflation_type[lx->depth] != 0){
    //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
    //}
    lx = lx->next_level;
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  h->hutch_compute_one_sample = g5_3D_hutchinson_mlmc_coarsest_PRECISION;
  
  if (g.probing) {
  for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++){
    if(g.my_rank==0)
      printf("\nEstimating trace at color %d\n", g.coloring_count);
      
    estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
    trace += estimate.acc_trace / estimate.sample_size;
    }
  } else {
      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
      trace += estimate.acc_trace / estimate.sample_size;
  }
  //if deflation vectors are available
  //if(g.trace_deflation_type[lx->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
  //}

  return trace;
}

complex_PRECISION g5_3D_connected_mlmc_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i, j;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);

  h->lx_i = h->finest_level;
  h->lx_j = h->finest_level;

  for(i = 0; i < g.num_levels; i++){
    h->lx_j = h->finest_level;
    for(j = 0; j < g.num_levels; j++){
    
      if(g.my_rank == 0) printf("\n Computing trace for C_{%d, %d} \n", i,j);
      
      h->hutch_compute_one_sample = connected_outer_PRECISION;
      
      // Call to blind from t+t’ with t’ = 0, ..., T-1
      for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected<g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {
        if (g.probing) {
          for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++) {
            if(g.my_rank==0) {
              printf("\nEstimating trace at color %d\n", g.coloring_count);
            }
            // TODO: should we indeed pass finest lvl? or some other? 
            estimate = hutchinson_blind_PRECISION(h->finest_level, h, 0, threading);
            trace += estimate.acc_trace / estimate.sample_size;
          }
        } else {
          estimate = hutchinson_blind_PRECISION(h->finest_level, h, 0, threading);
          trace += estimate.acc_trace / estimate.sample_size;
        }
      }
      
      
      int nlevs = g.num_levels;
      int idx = h->lx_i->depth * nlevs + h->lx_j->depth;  
      if (g.my_rank == 0) printf("\n Trace for C_{%d, %d} = %f\n, variance = %f", h->lx_i->depth, h->lx_j->depth, estimate.acc_trace / estimate.sample_size, g.variances[idx]);
      
      if (j < g.num_levels - 1)    
        h->lx_j = h->lx_j->next_level;
    }
      if (i < g.num_levels - 1)
        h->lx_i = h->lx_i->next_level;
  }

  
  
}


/* complex_PRECISION g5_3D_connected_mlmc_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;

  // for all but coarsest level
  lx = l;
  //for( i=0; i<1; i++ ) {
  for( i=0; i<g.num_levels-1; i++ ){
    // set the pointer to the mlmc difference operator
    h->hutch_compute_one_sample = g5_3D_connected_mlmc_difference_PRECISION;  
    
    for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected<g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {
      if (g.probing) {
        for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++) {
          if(g.my_rank==0) {
            printf("\nEstimating trace at color %d\n", g.coloring_count);
          }

          estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
          trace += estimate.acc_trace / estimate.sample_size;
        }
      } else {
        estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
        trace += estimate.acc_trace / estimate.sample_size;
      }
    }

    //if deflation vectors are available
    //if(g.trace_deflation_type[lx->depth] != 0){
    //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
    //}
    lx = lx->next_level;
  }

//  // coarsest level
//  // set the pointer to the coarsest-level Hutchinson estimator
//  h->hutch_compute_one_sample = g5_3D_hutchinson_mlmc_coarsest_PRECISION;

//  if (g.probing) {
//    for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++) {
//      if(g.my_rank==0) {
//        printf("\nEstimating trace at color %d\n", g.coloring_count);
//      }

//      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
//      trace += estimate.acc_trace / estimate.sample_size;
//    }
//  } else {
//      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
//      trace += estimate.acc_trace / estimate.sample_size;
//  }
//  //if deflation vectors are available
//  //if(g.trace_deflation_type[lx->depth] != 0){
//  //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
//  //}

  return trace;
}
*/

complex_PRECISION g5_3D_connected_split_orthogonal_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }

  // CONNECTED STEP #1. apply \Pi_{t+t'}
  // CONNECTED STEP #2. apply G5

  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
    //vector_PRECISION_ghg(  h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );
  }

  // CONNECTED STEP #3. apply level orthogonal difference operator
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Orthogonal difference (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
    {
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
      apply_solver_PRECISION( l, threading );
      
      // P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_R_PRECISION( h->mlmc_b1, p->x, l, threading );
      // P P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_P_PRECISION( h->mlmc_b2, h->mlmc_b1, l, threading );
      // (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
      vector_PRECISION_minus( h->mlmc_b1,  p->x, h->mlmc_b2, start, end, l );
    }
    // Prolongate solution to finest level
    //vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }
  }

  // CONNECTED STEP #4. apply \Pi_{t'}
  // CONNECTED STEP #5. apply G5
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // apply \Pi_{t'}
    //vector_PRECISION_ghg(  p->x, 0, l->inner_vector_size, l );
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    //g.time_slice = g.time_slice%g.global_lattice[0][0];
    //vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
    vector_PRECISION_ghg( h->mlmc_b1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->mlmc_b1, finest_l, threading );
  }

  // CONNECTED STEP #6. apply level orthogonal difference operator
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Orthogonal difference (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
    {
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
      apply_solver_PRECISION( l, threading );
      
      // P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_R_PRECISION( h->mlmc_b1, p->x, l, threading );
      // P P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_P_PRECISION( h->mlmc_b2, h->mlmc_b1, l, threading );
      // (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
      vector_PRECISION_minus( h->mlmc_b1,  p->x, h->mlmc_b2, start, end, l );
    }
    // Prolongate solution to finest level
    //vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }
  }
  
  
  // do the dot product AT FINEST
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);

    if ( type_appl==-1 ) {
      //if(g.trace_deflation_type[l->depth] != 0){
        //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
      //}
      aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, finest_l, threading );
    } else {
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
  }

  return aux;
}



complex_PRECISION g5_3D_connected_split_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx=0;

  // for all but coarsest level
  // lx = l;
  //for( i=0; i<g.num_levels-1 ;i++ ){
  //  if(g.my_rank == 0)
  //      printf("\nTrace at level %d split full rank operator\n", i+1);

    // set the pointer to the split full rank operator
    // h->hutch_compute_one_sample = hutchinson_split_intermediate_PRECISION;
    
    //if (g.probing) {
    //for (g.coloring_count = 1; g.coloring_count < g.num_colors[i+1] + 1; g.coloring_count++){
    //    estimate = hutchinson_blind_PRECISION(lx, h, 1, threading);
    //    trace += estimate.acc_trace / estimate.sample_size;
    // }
    //if(g.my_rank == 0){
    //    printf("\nTrace at level %d split full rank operator, Variance = %f\n", i+1, g.variances[i]);
    //    g.variances[i]=0.0;
    //}
    //} else {
    //    estimate = hutchinson_blind_PRECISION(lx, h, 1, threading);
    //    trace += estimate.acc_trace / estimate.sample_size;
    // }

    // if deflation vectors are available
  // if( g.trace_deflation_type[lx->depth] != 0 ){
  //  if( g.trace_deflation_type[lx->depth]==3 || g.trace_deflation_type[lx->depth]==5 ){
  //     trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
  //   }
  // }
  // lx = lx->next_level;    
// }

  // for all but coarsest level
  lx = l;
  //for( i=0; i<g.num_levels-1;i++ ){
  for( i=0; i < 1 ;i++ ){

    if(g.my_rank == 0){
      printf("\nTrace at level %d split orthogonal operator\n", i+1);
    }
    // set the pointer to the split orthogonal operator
    h->hutch_compute_one_sample = g5_3D_connected_split_orthogonal_PRECISION;
    
    for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected<g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {
      if (g.probing) {
      for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++){
          estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
          trace += estimate.acc_trace / estimate.sample_size;
      }
      if(g.my_rank == 0)
          printf("\nTrace at level %d split orthogonal operator, Variance = %f\n", i+1, g.variances[i]);
      } 
      else {
        estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
        trace += estimate.acc_trace / estimate.sample_size;
      }

    // if deflation vectors are available
    //if( g.trace_deflation_type[lx->depth] != 0 ){
    //  if( g.trace_deflation_type[lx->depth]==4 || g.trace_deflation_type[lx->depth]==5 ){
    //    trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
    //  }
    // }
    }
      lx = lx->next_level; 
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  //if(g.my_rank == 0)
  //       printf("\nTrace at level %d \n", g.num_levels);

  //h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  //if (g.probing) {
  //  for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++){
  //      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
  //      trace += estimate.acc_trace / estimate.sample_size;
  //   }
  //   if(g.my_rank == 0)
  //      printf("\nTrace at level %d split orthogonal operator, Variance = %f\n", g.num_levels, g.variances[g.num_levels-1]);
  //  } else {
  //      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
  //      trace += estimate.acc_trace / estimate.sample_size;
  //  }

  // if deflation vectors are available
  //if(g.trace_deflation_type[lx->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
  // }
  return trace;
}








void connected_mlmc_PRECISION_non_difference( vector_PRECISION out, vector_PRECISION in, level_struct *l,  hutchinson_PRECISION_struct* h, struct Thread *threading ){

  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }
  
  // Apply Gamma_5 at the finest
  {    
    gamma5_PRECISION( h->mlmc_testing, in, finest_l, threading );
  }
  
  // Solve at current level : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    //if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    //} else {
    // // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    // }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level
    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(out, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  

}







void connected_mlmc_PRECISION_difference( vector_PRECISION out, vector_PRECISION in, level_struct *l,  hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }
  
  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );
    
  // Apply Gamma_5 at the finest
  {    
    gamma5_PRECISION( h->mlmc_testing, in, finest_l, threading );
  }
  
  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    //if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    //} else {
    // // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    // }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level
    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  
  // LEVEL DIFFERENCE : SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict (iteratively from finest)
  // 2. invert
  // 3. Prolongate (iteratively to finest)
  {
    int start, end, start_tmp, end_tmp, start_tmp_c, end_tmp_c;
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    //if ( type_appl==-1 ) {
      // Copy gamma_5 x at the finest
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

      // Restrict it from finest to NEXT level
      for (int i = 0; i < l_index+1; ++i) {
        level_struct *fine   = levels[i];
        level_struct *coarse = levels[i + 1];

        compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
        compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
        // restrict from fine to coarse
        apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
        vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

        if (g.my_rank == 0)
          printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
      }

      // Copy restricted vector to rhs at the next level
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
      vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, end_tmp, l->next_level);
  // } else {
    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
  // }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );

    // Prolongate Solution from NEXT level to finest
    vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
    for (int i = l_index+1; i > 0; --i) {
      level_struct *coarse = levels[i];
      // finer level (towards finest)
      level_struct *fine   = levels[i - 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // prolongate from coarse to fine
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      // copy result to fine level
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

      if(g.my_rank == 0) {
        printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in second term vector
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    //copy for dot product
    vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  
  // subtract the results
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    vector_PRECISION_minus( out, first_term, second_term, start, end, finest_l);
  }
  
  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );  

}



//TODO: better naming (?)
// builds C_ij by calling mlmc_difference and mlmc_non_difference
complex_PRECISION connected_outer_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

  complex_PRECISION trace = 0.0;
  struct sample estimate;
  
  level_struct *finest_l = h->finest_level;
  level_struct* lx_i;
  level_struct* lx_j;
  lx_i = h->lx_i;
  lx_j = h->lx_j;
  
  int i = lx_i->depth;
  int j = lx_j->depth;
  
  // Allocate buffer at finest
  vector_PRECISION op1 = NULL;
  PUBLIC_MALLOC( op1, complex_PRECISION, finest_l->inner_vector_size );
  
  // for all but coarsest level
  //lx = l;
  
  // CONNECTED STEP #1. apply \Pi_{t+t'}   (\Pi_{t+t'} x)
  {
    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
  }

  // CONNECTED STEP #2. apply Operator S_j    (S_j \Pi_{t+t'} x)
  {
    if( j == g.num_levels - 1){   
      connected_mlmc_PRECISION_non_difference( op1, h->rademacher_vector, lx_j,  h, threading );
    }
    else{
      connected_mlmc_PRECISION_difference( op1, h->rademacher_vector,  lx_j,  h, threading );
    }
  }
  
  // CONNECTED STEP #3. apply \Pi{t’}    (\Pi_{t'} S_j \Pi_{t+t'} x) 
  {
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    vector_PRECISION_ghg( op1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
  }
  
  // CONNECTED STEP #4. apply Difference operator S_i    (S_i \Pi_{t'} S_j \Pi_{t+t'} x)
  {
    if( i == g.num_levels - 1){   
      connected_mlmc_PRECISION_non_difference( op1, op1,  lx_i,  h, threading );
    }
    else{
      connected_mlmc_PRECISION_difference( op1, op1,  lx_i,  h, threading );
    }
  }

  
  // CONNECTED STEP #5 do the dot product AT FINEST       (x^H \Pi_{t+t'} S_i \Pi_{t'} S_j \Pi_{t+t'} x))
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);

    //if ( type_appl==-1 ) {
    //  if(g.trace_deflation_type[l->depth] != 0){
    //    hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
    //  }
    trace = global_inner_product_PRECISION( h->rademacher_vector, op1, p->v_start, p->v_end, finest_l, threading );
    //  } else {
    //  aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    // }
  }
  
  PUBLIC_FREE( op1, complex_PRECISION, finest_l->inner_vector_size );
  
  return trace;

}









/*

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

*/
