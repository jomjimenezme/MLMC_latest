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


void cpu_fgmres_PRECISION_struct_init( gmres_PRECISION_struct *p ) {

/*********************************************************************************
* Initializes all declared pointers with NULL.                              
*********************************************************************************/

  p->Z = NULL;
  p->V = NULL;
  p->H = NULL;
  p->x = NULL;
  p->b = NULL;
  p->r = NULL;
  p->w = NULL;
  p->y = NULL;
  p->gamma = NULL;
  p->c = NULL;
  p->s = NULL;
  p->shift = 0;
  p->preconditioner = NULL;
  p->eval_operator = NULL;
#ifdef GCR_SMOOTHER
  p->use_gcr = 0;
#endif
#ifdef RICHARDSON_SMOOTHER
  p->use_richardson = 0;
  p->richardson_update_omega = 1;
#endif

  // copy of Hesselnberg matrix
#if defined(GCRODR) && defined(POLYPREC)
  p->gcrodr_PRECISION.eigslvr.Hc = NULL;
  p->polyprec_PRECISION.eigslvr.Hc = NULL;
#elif defined(GCRODR)
  p->gcrodr_PRECISION.eigslvr.Hc = NULL;
#elif defined(POLYPREC)
  p->polyprec_PRECISION.eigslvr.Hc = NULL;
#endif

#ifdef POLYPREC
  p->polyprec_PRECISION.Hcc = NULL; 
  p->polyprec_PRECISION.L = NULL;
  p->polyprec_PRECISION.col_prods = NULL;
  p->polyprec_PRECISION.accum_prod = NULL;
  p->polyprec_PRECISION.product = NULL;
  p->polyprec_PRECISION.temp = NULL;
  p->polyprec_PRECISION.h_ritz = NULL;
  p->polyprec_PRECISION.lejas = NULL;
  p->polyprec_PRECISION.random_rhs = NULL;
  p->polyprec_PRECISION.xtmp = NULL;

  p->polyprec_PRECISION.eigslvr.vl = NULL;
  p->polyprec_PRECISION.eigslvr.vr = NULL;
  p->polyprec_PRECISION.dirctslvr.ipiv = NULL;
  p->polyprec_PRECISION.dirctslvr.x = NULL;
  p->polyprec_PRECISION.dirctslvr.b = NULL;
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->Va = NULL;
  p->Za = NULL;
#endif

//#ifdef BLOCK_JACOBI
#if 0
  p->block_jacobi_PRECISION.b_backup = NULL;
  p->block_jacobi_PRECISION.xtmp = NULL;
  local_fgmres_PRECISION_struct_init( &(p->block_jacobi_PRECISION.local_p) );
#endif

#ifdef RICHARDSON_SMOOTHER
  p->omega = NULL;
#endif

#ifdef CUDA_OPT
  p->w_componentwise_gpu = NULL;
  p->r_componentwise_gpu = NULL;
  p->x_componentwise_gpu = NULL;
  p->x_gpu = NULL;
  p->b_componentwise_gpu = NULL;
  p->b_gpu = NULL;
#endif

#ifdef GCR_SMOOTHER
  p->gcr_buffer_dotprods = NULL;
  p->gcr_betas_dotprods  = NULL;
#endif

  p->print_iters = 0;
}


void cpu_fgmres_PRECISION_struct_alloc( int m, int n, int vl, PRECISION tol, const int type, const int prec_kind,
                                    void (*precond)(), void (*eval_op)(), gmres_PRECISION_struct *p, level_struct *l ) {

/*********************************************************************************
* Allocates memory for the fgmres struct and sets its values.                  
* int m: Restart length                                                            
* int n: Number of restarts                                                        
* int vl: System size                                                              
* PRECISION tol: Tolerance for relative residual                                         
* const int type: Specifies the problem for which fgmres should be applied               
*                 (_GLOBAL_FGMRES, _K_CYCLE, _COARSE_GMRES)                              
* const int prec_kind: type of preconditioning: _RIGHT (flexible preconditioner),
*                                               _LEFT (stationary preconditioner)
*                                               or _NOTHING                                     
* void (*precond): Function pointer to the preconditioner                           
*********************************************************************************/  
  
  long int total=0; 
  int i, k=0;
  
  p->restart_length = m;
  p->num_restart = n;
  p->preconditioner = precond;
  p->eval_operator = eval_op; 
  p->tol = tol;
  p->kind = prec_kind;
  
  if(m > 0) {
  total += (m+1)*m; // Hessenberg matrix
  MALLOC( p->H, complex_PRECISION*, m );
  
  total += (5+m)*vl; // x, r, b, w, V
  MALLOC( p->V, complex_PRECISION*, m+1 );
  
  if ( precond != NULL ) {
    if ( prec_kind == _RIGHT ) {
      total += (m+1)*vl; // Z
      k = m+1;
    } else {
      total += vl;
      k = 1;
    }
    MALLOC( p->Z, complex_PRECISION*, k );
  } else {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      total += (m+2)*vl;
      k = m+2;
      MALLOC( p->Z, complex_PRECISION*, k );
    }
#else
#ifdef GCR_SMOOTHER
    if ( l->depth==0 ) {
      total += (m+1)*vl; // Z
      k = m+1;
      MALLOC( p->Z, complex_PRECISION*, k );
    }
    else { k=0; }
#else
    k = 0;
#endif
#endif
  }
  
  total += 4*(m+1); // y, gamma, c, s
  
  p->H[0] = NULL; // allocate connected memory
  MALLOC( p->H[0], complex_PRECISION, total );
  
  p->total_storage = total;
  total = 0;
  
  // ordering: H, y, gamma, c, s, w, V, Z, x, r, b
  // H
  for ( i=1; i<m; i++ )
    p->H[i] = p->H[0] + i*(m+1);
  total += m*(m+1);
  
  // y
  p->y = p->H[0] + total; total += m+1;
  // gamma
  p->gamma = p->H[0] + total; total += m+1;
  // c
  p->c = p->H[0] + total; total += m+1;
  // s
  p->s = p->H[0] + total; total += m+1;
  // w
  p->w = p->H[0] + total; total += vl;
  // V
  for ( i=0; i<m+1; i++ ) {
    p->V[i] = p->H[0] + total; total += vl;
  }
  // Z
  for ( i=0; i<k; i++ ) {
    p->Z[i] = p->H[0] + total; total += vl;
  }

  // x
  p->x = p->H[0] + total; total += vl;
  // r
  p->r = p->H[0] + total; total += vl;
  // b
  p->b = p->H[0] + total; total += vl;
  
  ASSERT( p->total_storage == total );
  }
  
  if ( type == _GLOBAL_FGMRES ) {    
    p->timing = 1;
    p->print = g.vt.evaluation?0:1;
    p->initial_guess_zero = 1;
    p->shift = 0;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(g.op_PRECISION);
  } else if ( type == _K_CYCLE ) {
    // these settings also work for GMRES as a smoother
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->shift = 0;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(l->s_PRECISION.op);
  } else if ( type == _COARSE_GMRES ) {
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->layout = -1;
    p->shift = 0;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    if ( g.odd_even )
      p->op = &(l->oe_op_PRECISION);
    else  
      p->op = &(l->s_PRECISION.op);
  } else {
    ASSERT( type < 3 );
  }

#if defined(GCRODR) || defined(POLYPREC)
  if (l->level==0) {
#endif

  // FIXME : is this function-pointer-assignment really necessary ?
#if defined(GCRODR) || defined(POLYPREC)
  //p->polyprec_PRECISION.eigslvr.eigslvr_PRECISION = eigslvr_PRECISION;
#endif

  // copy of Hesselnberg matrix
#if defined(GCRODR) && defined(POLYPREC)
  MALLOC(p->gcrodr_PRECISION.eigslvr.Hc, complex_PRECISION*, m);
  p->polyprec_PRECISION.eigslvr.Hc = p->gcrodr_PRECISION.eigslvr.Hc;
  p->gcrodr_PRECISION.eigslvr.Hc[0] = NULL; // allocate connected memory
  MALLOC( p->gcrodr_PRECISION.eigslvr.Hc[0], complex_PRECISION, m*(m+1) );
  for ( i=1; i<m; i++ )
    p->gcrodr_PRECISION.eigslvr.Hc[i] = p->gcrodr_PRECISION.eigslvr.Hc[0] + i*(m+1);
  p->polyprec_PRECISION.eigslvr.Hc[0] = p->gcrodr_PRECISION.eigslvr.Hc[0];
#elif defined(GCRODR)
  MALLOC(p->gcrodr_PRECISION.eigslvr.Hc, complex_PRECISION*, m);
  p->gcrodr_PRECISION.eigslvr.Hc[0] = NULL; // allocate connected memory
  MALLOC( p->gcrodr_PRECISION.eigslvr.Hc[0], complex_PRECISION, m*(m+1) );
  for ( i=1; i<m; i++ )
    p->gcrodr_PRECISION.eigslvr.Hc[i] = p->gcrodr_PRECISION.eigslvr.Hc[0] + i*(m+1);
#elif defined(POLYPREC)
  MALLOC(p->polyprec_PRECISION.eigslvr.Hc, complex_PRECISION*, m);
  p->polyprec_PRECISION.eigslvr.Hc[0] = NULL; // allocate connected memory
  MALLOC( p->polyprec_PRECISION.eigslvr.Hc[0], complex_PRECISION, m*(m+1) );
  for ( i=1; i<m; i++ )
    p->polyprec_PRECISION.eigslvr.Hc[i] = p->polyprec_PRECISION.eigslvr.Hc[0] + i*(m+1);
#endif

#ifdef POLYPREC
  int d_max = (g.polyprec_d_setup>g.polyprec_d_solve)?g.polyprec_d_setup:g.polyprec_d_solve;
  p->polyprec_PRECISION.d_poly = d_max;
  int d_poly = p->polyprec_PRECISION.d_poly;

  MALLOC( p->polyprec_PRECISION.col_prods, complex_PRECISION, d_poly);
  MALLOC( p->polyprec_PRECISION.h_ritz, complex_PRECISION, d_poly);
  MALLOC( p->polyprec_PRECISION.lejas, complex_PRECISION, d_poly);
  MALLOC( p->polyprec_PRECISION.random_rhs, complex_PRECISION, vl );
  MALLOC( p->polyprec_PRECISION.accum_prod, complex_PRECISION, vl );
  MALLOC( p->polyprec_PRECISION.product, complex_PRECISION, vl );
  MALLOC( p->polyprec_PRECISION.temp, complex_PRECISION, vl );

  MALLOC( p->polyprec_PRECISION.xtmp, complex_PRECISION, vl );

  MALLOC( p->polyprec_PRECISION.Hcc, complex_PRECISION, d_poly*d_poly );
  MALLOC( p->polyprec_PRECISION.L, complex_PRECISION*, d_poly+ 1);

  p->polyprec_PRECISION.L[0] = NULL;

  MALLOC( p->polyprec_PRECISION.L[0], complex_PRECISION, (d_poly+1)*d_poly );

  for (i=1; i<d_poly+1; i++)
  {
    p->polyprec_PRECISION.L[i] = p->polyprec_PRECISION.L[0] + i*d_poly;
  }

  MALLOC( p->polyprec_PRECISION.dirctslvr.ipiv, int, d_poly);
  MALLOC( p->polyprec_PRECISION.dirctslvr.x, complex_PRECISION, d_poly);
  MALLOC( p->polyprec_PRECISION.dirctslvr.b, complex_PRECISION, d_poly);

  p->polyprec_PRECISION.dirctslvr.N = d_poly;
  p->polyprec_PRECISION.dirctslvr.lda = d_poly; // m here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  p->polyprec_PRECISION.dirctslvr.ldb = d_poly;
  p->polyprec_PRECISION.dirctslvr.nrhs = 1;
  p->polyprec_PRECISION.dirctslvr.Hcc = p->polyprec_PRECISION.Hcc;
  p->polyprec_PRECISION.dirctslvr.dirctslvr_PRECISION = dirctslvr_PRECISION;

  MALLOC( p->polyprec_PRECISION.eigslvr.vl, complex_PRECISION, d_poly*d_poly );
  MALLOC( p->polyprec_PRECISION.eigslvr.vr, complex_PRECISION, d_poly*d_poly );

  p->polyprec_PRECISION.eigslvr.jobvl = 'N';
  p->polyprec_PRECISION.eigslvr.jobvr = 'N';

  p->polyprec_PRECISION.eigslvr.N = d_poly;
  p->polyprec_PRECISION.eigslvr.lda = p->restart_length + 1;
  p->polyprec_PRECISION.eigslvr.ldvl = d_poly;
  p->polyprec_PRECISION.eigslvr.ldvr = d_poly;
  p->polyprec_PRECISION.eigslvr.w = p->polyprec_PRECISION.h_ritz;
  p->polyprec_PRECISION.Hc = p->polyprec_PRECISION.eigslvr.Hc;
  p->polyprec_PRECISION.eigslvr.eigslvr_PRECISION = eigslvr_PRECISION;    

  p->polyprec_PRECISION.update_lejas = 1;
  p->polyprec_PRECISION.preconditioner = NULL;
  p->polyprec_PRECISION.preconditioner_bare = p->preconditioner;
  p->polyprec_PRECISION.syst_size = vl;

  p->polyprec_PRECISION.eigslvr.A = p->polyprec_PRECISION.Hc[0];
#endif

#if defined(GCRODR) || defined(POLYPREC)
  }
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->syst_size = vl;
  MALLOC( p->Va, complex_PRECISION*, m+2 );
  MALLOC( p->Za, complex_PRECISION*, m+2 );
  p->Va[0] = NULL;
  p->Za[0] = NULL;
  MALLOC( p->Va[0], complex_PRECISION, (m+2)*vl );
  MALLOC( p->Za[0], complex_PRECISION, (m+2)*vl );

  for ( i=1; i<m+2; i++ )
  {
    p->Va[i] = p->Va[0] + i*vl;
    p->Za[i] = p->Za[0] + i*vl;
  }

#ifdef PERS_COMMS
  g.pers_comms_nrZas = m+2;
#endif
#endif

//#ifdef BLOCK_JACOBI
#if 0
  p->block_jacobi_PRECISION.syst_size = vl;

  if (l->level==0) {
    // these two always go together
    p->block_jacobi_PRECISION.BJ_usable = 0;
    p->block_jacobi_PRECISION.local_p.polyprec_PRECISION.update_lejas = 1;

    MALLOC( p->block_jacobi_PRECISION.b_backup, complex_PRECISION, vl );
    MALLOC( p->block_jacobi_PRECISION.xtmp, complex_PRECISION, vl );

    p->block_jacobi_PRECISION.local_p.polyprec_PRECISION.d_poly = g.local_polyprec_d;

    local_fgmres_PRECISION_struct_alloc( g.local_polyprec_d, 1, vl, g.coarse_tol, 
                                         _COARSE_GMRES, _NOTHING, NULL,
                                         coarse_local_apply_schur_complement_PRECISION,
                                         &(p->block_jacobi_PRECISION.local_p), l );
  }
#endif

#ifdef RICHARDSON_SMOOTHER
  p->richardson_sub_degree = g.richardson_sub_degree;
  MALLOC( p->omega, PRECISION, p->richardson_sub_degree );
  p->richardson_factor = g.richardson_factor;
#endif

#ifdef CUDA_OPT
  // if these conditions are fulfilled, this is then the finest-level smoother
  // with method=4, which is GCR, GMRES or Richardson
  if ( g.method==4 && prec_kind==_NOTHING && l->depth==0 ) {
    p->gpu_syst_size = vl;
    CUDA_MALLOC( p->w_componentwise_gpu, cu_cmplx_PRECISION, vl );
    CUDA_MALLOC( p->r_componentwise_gpu, cu_cmplx_PRECISION, vl );
    CUDA_MALLOC( p->x_componentwise_gpu, cu_cmplx_PRECISION, vl );
    CUDA_MALLOC( p->x_gpu, cu_cmplx_PRECISION, vl );
    CUDA_MALLOC( p->b_componentwise_gpu, cu_cmplx_PRECISION, vl );
    CUDA_MALLOC( p->b_gpu, cu_cmplx_PRECISION, vl );
  }
#endif

#ifdef GCR_SMOOTHER
  MALLOC( p->gcr_buffer_dotprods, complex_PRECISION, p->restart_length+1 );
  MALLOC( p->gcr_betas_dotprods,  complex_PRECISION, p->restart_length+1 );
#endif
}


void cpu_fgmres_PRECISION_struct_free( gmres_PRECISION_struct *p, level_struct *l ) {

/*********************************************************************************
* Frees the allocated space for the gmres struct p.                            
*********************************************************************************/   
  int k=0;
  
  if ( p->preconditioner != NULL ) {
    if ( p->kind == _RIGHT )
      k = p->restart_length+1;
    else
      k = 1;
  } else {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      k = p->restart_length+2;
    }
#else
#ifdef GCR_SMOOTHER
    if ( l->depth==0 ) {
      k = p->restart_length+1;
    }
    else { k=0; }
#else
    k = 0;
#endif
#endif
  }

  if(p->restart_length > 0) {
  FREE( p->H[0], complex_PRECISION, p->total_storage );
  FREE( p->H, complex_PRECISION*, p->restart_length );
  FREE( p->V, complex_PRECISION*, p->restart_length+1 );
  
  if ( p->Z != NULL )
    FREE( p->Z, complex_PRECISION*, k );
  }
  
  p->D = NULL;
  p->clover = NULL;

  // --- COARSEST-LEVEL IMPROVEMENTS

#if defined(GCRODR) || defined(POLYPREC)
  if (l->level==0) {
#endif

  // copy of Hesselnberg matrix
#if defined(GCRODR) && defined(POLYPREC)
  int m = p->restart_length;
  FREE( p->gcrodr_PRECISION.eigslvr.Hc[0], complex_PRECISION, m*(m+1) );
  FREE(p->gcrodr_PRECISION.eigslvr.Hc, complex_PRECISION*, m);
#elif defined(GCRODR)
  int m = p->restart_length;
  FREE( p->gcrodr_PRECISION.eigslvr.Hc[0], complex_PRECISION, m*(m+1) );
  FREE(p->gcrodr_PRECISION.eigslvr.Hc, complex_PRECISION*, m);
#elif defined(POLYPREC)
  int m = p->restart_length;
  FREE( p->polyprec_PRECISION.eigslvr.Hc[0], complex_PRECISION, m*(m+1) );
  FREE(p->polyprec_PRECISION.eigslvr.Hc, complex_PRECISION*, m);
#endif

#ifdef POLYPREC
  int d_poly = p->polyprec_PRECISION.d_poly;
  int vl = p->polyprec_PRECISION.syst_size;
  FREE( p->polyprec_PRECISION.Hcc, complex_PRECISION, d_poly*d_poly );
  FREE( p->polyprec_PRECISION.L[0], complex_PRECISION, (d_poly+1)*d_poly );
  FREE( p->polyprec_PRECISION.L, complex_PRECISION*, d_poly+1 );
  FREE( p->polyprec_PRECISION.h_ritz,complex_PRECISION, d_poly );
  FREE( p->polyprec_PRECISION.lejas,complex_PRECISION, d_poly );
  FREE( p->polyprec_PRECISION.accum_prod, complex_PRECISION, vl );
  FREE( p->polyprec_PRECISION.product, complex_PRECISION, vl );    
  FREE( p->polyprec_PRECISION.temp, complex_PRECISION, vl );    
  FREE( p->polyprec_PRECISION.xtmp, complex_PRECISION, vl );
  FREE( p->polyprec_PRECISION.random_rhs, complex_PRECISION, vl );
  FREE( p->polyprec_PRECISION.col_prods, complex_PRECISION, d_poly );

  FREE( p->polyprec_PRECISION.eigslvr.vl,complex_PRECISION, d_poly*d_poly );
  FREE( p->polyprec_PRECISION.eigslvr.vr,complex_PRECISION, d_poly*d_poly );  

  FREE( p->polyprec_PRECISION.dirctslvr.ipiv, int, d_poly );
  FREE( p->polyprec_PRECISION.dirctslvr.x, complex_PRECISION, d_poly );
  FREE( p->polyprec_PRECISION.dirctslvr.b, complex_PRECISION, d_poly );
#endif

#if defined(GCRODR) || defined(POLYPREC)
  }
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    FREE( p->Va[0], complex_PRECISION*, (p->restart_length+2)*p->syst_size );
    FREE( p->Za[0], complex_PRECISION*, (p->restart_length+2)*p->syst_size );
    FREE( p->Va, complex_PRECISION, p->restart_length+2 );
    FREE( p->Za, complex_PRECISION, p->restart_length+2 );
#endif

//#ifdef BLOCK_JACOBI
#if 0
  if (l->level==0) {
    FREE( p->block_jacobi_PRECISION.b_backup, complex_PRECISION, p->block_jacobi_PRECISION.syst_size );
    FREE( p->block_jacobi_PRECISION.xtmp, complex_PRECISION, p->block_jacobi_PRECISION.syst_size );

    local_fgmres_PRECISION_struct_free( &(p->block_jacobi_PRECISION.local_p), l );
  }
#endif

#ifdef RICHARDSON_SMOOTHER
  FREE( p->omega, PRECISION, p->richardson_sub_degree );
#endif

#ifdef CUDA_OPT
  // if these conditions are fulfilled, this is then the finest-level smoother
  // with method=4, which is GCR, GMRES or Richardson
  if ( g.method==4 && p->kind==_NOTHING && l->depth==0 ) {
    CUDA_FREE( p->w_componentwise_gpu, cu_cmplx_PRECISION, p->gpu_syst_size );
    CUDA_FREE( p->r_componentwise_gpu, cu_cmplx_PRECISION, p->gpu_syst_size );
    CUDA_FREE( p->x_componentwise_gpu, cu_cmplx_PRECISION, p->gpu_syst_size );
    CUDA_FREE( p->x_gpu, cu_cmplx_PRECISION, p->gpu_syst_size );
    CUDA_FREE( p->b_componentwise_gpu, cu_cmplx_PRECISION, p->gpu_syst_size );
    CUDA_FREE( p->b_gpu, cu_cmplx_PRECISION, p->gpu_syst_size );
  }
#endif

#ifdef GCR_SMOOTHER
  FREE( p->gcr_buffer_dotprods, complex_PRECISION, p->restart_length+1 );
  FREE( p->gcr_betas_dotprods,  complex_PRECISION, p->restart_length+1 );
#endif
}


int fgmres_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Uses FGMRES to solve the system D x = b, where b is taken from p->b and x is 
* stored in p->x.                                                              
*********************************************************************************/  
  RangeHandleType profilingRangeResFgmres = startProfilingRange("Restarted FGMRES (PRECISION)");
  // RE-ENABLE !
  //printf0("WITHIN fgmres_PRECISION(...), depth=%d \n", l->depth);

  // start and end indices for vector functions depending on thread
  int start;
  int end;

  int j=-1, finish=0, iter=0, il, ol, res;
  complex_PRECISION gamma0 = 0;

  complex_PRECISION beta = 0;

  double norm_r0=1, gamma_jp1=1, t0=0, t1=0;
  START_LOCKED_MASTER(threading)

  if ( l->depth==0 && ( p->timing || p->print ) ) prof_init( l );

  if ( l->level==0 && g.num_levels > 1 && g.interpolation ) p->tol = g.coarse_tol;
  if ( l->depth > 0 ) p->timing = 1;
  if ( l->depth == 0 ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( p->print && g.print > 0 ) printf0("+----------------------------------------------------------+\n");
#endif
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  // Creation of CUDA Events, for time measurement and GPU sync
  //cudaEvent_t start_event_copy, stop_event_copy, start_event_comp, stop_event_comp;
  //cuda_safe_call( cudaEventCreate(&start_event_copy) );
  //cuda_safe_call( cudaEventCreate(&stop_event_copy) );
  //cuda_safe_call( cudaEventCreate(&start_event_comp) );
  //cuda_safe_call( cudaEventCreate(&stop_event_comp) );
  //cudaStream_t *streams_gmres = p->streams;

  // RE-ENABLE !
  //printf0("p->num_restart = %d\n", p->num_restart);

  for( ol=0; ol<p->num_restart && finish==0; ol++ )  {
    RangeHandleType profilingRangeFgmres = startProfilingRange("FGMRES (PRECISION)");

    // RE-ENABLE !
    ///printf0("for loop of fgmres_PRECISION, iter=%d, depth=%d \n", ol, l->depth);
  
    if( ol == 0 && p->initial_guess_zero ) {
      res = _NO_RES;
      vector_PRECISION_copy( p->r, p->b, start, end, l );
    } else {
      res = _RES;
      if ( p->kind == _LEFT && p->preconditioner ) {
        apply_operator_PRECISION( p->Z[0], p->x, p, l, threading );

        if ( p->shift ) vector_PRECISION_saxpy( p->Z[0], p->Z[0], p->x, p->shift, start, end, l );
        if ( g.method == 5 ) {
          START_LOCKED_MASTER(threading)
          g.bicgstab_tol = (!g.mixed_precision)?p->tol:MAX( 1E-3, (p->tol/(gamma_jp1/norm_r0))*5E-1 );
          END_LOCKED_MASTER(threading)
        }
        p->preconditioner( p->w, NULL, p->Z[0], _NO_RES, l, threading );
      } else {
        apply_operator_PRECISION( p->w, p->x, p, l, threading ); // compute w = D*x
      }
      vector_PRECISION_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
    }

    gamma0 = (complex_PRECISION) global_norm_PRECISION( p->r, p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)
    START_MASTER(threading)
    p->gamma[0] = gamma0;
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if( ol == 0) {
      norm_r0 = creal(p->gamma[0]);
    }
    
    vector_PRECISION_real_scale( p->V[0], p->r, 1/p->gamma[0], start, end, l ); // v_0 = r / gamma_0
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, 0, p->preconditioner, p->shift, p, l, threading );
    }
#endif   
    
    for( il=0; il<p->restart_length && finish==0; il++) {

      // RE-ENABLE !!
      //printf0("\ninner for loop of fgmres_PRECISION, iter=%d, depth=%d \n\n", il, l->depth);

      j = il; iter++;
      if ( g.method == 5 ) {
        START_LOCKED_MASTER(threading)
        g.bicgstab_tol = (!g.mixed_precision)?p->tol:MAX( 1E-3, (p->tol/(gamma_jp1/norm_r0))*5E-1 );
        END_LOCKED_MASTER(threading)
      }
      
      // one step of Arnoldi
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      if ( l->level == 0 && l->depth > 0 ) {
        if ( !arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j+1, p->preconditioner, p->shift, p, l, threading ) ) {
          printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+2, j+1 );
          break;
        }
      } else {
        if ( !arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p->shift, p, l, threading ) ) {
          printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
          break;
        }
      }
#else
      if ( !arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p->shift, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
        break;
      }
#endif
      
      if ( cabs( p->H[j][j+1] ) > p->tol/10 ) {
        qr_update_PRECISION( p->H, p->s, p->c, p->gamma, j, l, threading );
        gamma_jp1 = cabs( p->gamma[j+1] );
        
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
        if ( iter%10 == 0 || p->preconditioner != NULL || l->depth > 0 ) {
          START_MASTER(threading)
          if ( p->print && g.print > 0 )
            printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter, gamma_jp1/norm_r0 );
          END_MASTER(threading)
        }
#endif

        //printf0( "%.14f\n",gamma_jp1/norm_r0 );

        if( gamma_jp1/norm_r0 < p->tol || gamma_jp1/norm_r0 > 1E+5 ) { // if satisfied ... stop
//#ifdef BLOCK_JACOBI
#if 0
          if ( l->level==0 )
          {
            // backup of p->x, just in case tol hasn't been reached we need to restore ...
            vector_PRECISION_copy( p->block_jacobi_PRECISION.xtmp, p->x, start, end, l );

            compute_solution_PRECISION( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                        p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );

            p->eval_operator( p->w, p->x, p->op, l, threading );
            vector_PRECISION_minus( p->r, p->block_jacobi_PRECISION.b_backup, p->w, start, end, l ); // compute r = b - w
            PRECISION norm_r0xx = global_norm_PRECISION( p->block_jacobi_PRECISION.b_backup, start, end, l, threading );
            PRECISION betaxx = global_norm_PRECISION( p->r, start, end, l, threading );
            if ( betaxx/norm_r0xx < p->tol ) {
              finish = 1;
            } else {
              // restore p->x
              vector_PRECISION_copy( p->x, p->block_jacobi_PRECISION.xtmp, start, end, l );
            }
            START_MASTER(threading)
            if ( betaxx/norm_r0xx > 1E+5 ) printf0("Divergence of fgmres_PRECISION, iter = %d, level=%d\n", iter, l->level );
            END_MASTER(threading)
          } else {
            finish = 1;
            START_MASTER(threading)
            if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_PRECISION, iter = %d, level=%d\n", iter, l->level );
            END_MASTER(threading)
          }
#else
          finish = 1;
          START_MASTER(threading)
          if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_PRECISION, iter = %d, level=%d\n", iter, l->level );
          END_MASTER(threading)
#endif
        }
      } else {
        printf0("depth: %d, iter: %d, p->H(%d,%d) = %+lf+%lfi\n", l->depth, iter, j+1, j, CSPLIT( p->H[j][j+1] ) );
        finish = 1;
        break;
      }
    } // end of a single restart
//#ifdef BLOCK_JACOBI
#if 0
    if ( l->level==0 ) {
      if ( finish==0 ) {
        compute_solution_PRECISION( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                    p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );
      }
    } else {
      compute_solution_PRECISION( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                  p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );
    }
#else
    compute_solution_PRECISION( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );
#endif
    endProfilingRange(profilingRangeFgmres);
  } // end of fgmres
  
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.total_time = t1-t0; g.iter_count = iter; g.norm_res = gamma_jp1/norm_r0; }
  END_LOCKED_MASTER(threading)
  
  if ( p->print ) {
#ifdef FGMRES_RESTEST
    apply_operator_PRECISION( p->w, p->x, p, l, threading );
    vector_PRECISION_minus( p->r, p->b, p->w, start, end, l );
    beta = global_norm_PRECISION( p->r, p->v_start, p->v_end, l, threading );
#else
    beta = gamma_jp1;
#endif
    START_MASTER(threading)
    g.norm_res = creal(beta)/norm_r0;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if ( g.print > 0 ) printf0("+----------------------------------------------------------+\n\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|       FGMRES iterations: %-6d coarse average: %-6.2lf   |\n", iter,
            ((double)g.coarse_iter_count)/((double)iter) );
    printf0("| exact relative residual: ||r||/||b|| = %e      |\n", creal(beta)/norm_r0 );
    printf0("| elapsed wall clock time: %-8.4lf seconds                |\n", t1-t0 );
    if ( g.coarse_time > 0 ) 
      printf0("|        coarse grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.coarse_time, 100*(g.coarse_time/(t1-t0)) );
    printf0("|        coarsest grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.coarsest_time, 100*(g.coarsest_time/(t1-t0)) );
    printf0("| coarsest grid matmul time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.matmul_time, 100*(g.matmul_time/(t1-t0)) );
//#ifdef BLOCK_JACOBI
#if 0
    printf0("|     coarsest grid BJ time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.bj_time, 100*(g.bj_time/(t1-t0)) );
#endif
#ifdef GCRODR
    printf0("| coarsest grid GCRODR LSP time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.gcrodr_LSP_time, 100*(g.gcrodr_LSP_time/(t1-t0)) );
    printf0("| coarsest grid GCRODR AB time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.gcrodr_buildAB_time, 100*(g.gcrodr_buildAB_time/(t1-t0)) );
    printf0("| coarsest grid GCRODR CU time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.gcrodr_buildCU_time, 100*(g.gcrodr_buildCU_time/(t1-t0)) );
#endif
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
#ifdef CUDA_OPT
    printf0("|        max used mem/GPU: %-8.2le GB                     |\n", g.max_gpu_storage/1024.0 );
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }
  
#ifdef COARSE_RES
  if ( l->depth > 0 ) {
    START_MASTER(threading)
    char number[3]; sprintf( number, "%2d", 31+l->depth ); printf0("\033[1;%2sm|", number );
    printf0(" - depth: %d, gmres iter: %2d, approx rel res: %le |", l->depth, iter, gamma_jp1/norm_r0 );
    printf0("\033[0m\n"); fflush(0);
    END_MASTER(threading)
  }
#endif

  if ( l->level == 0 ) {
    START_LOCKED_MASTER(threading)
    g.coarse_iter_count += iter;
    END_LOCKED_MASTER(threading)
  }
    
  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      START_LOCKED_MASTER(threading)
      printf0("solve iter: %d\n", iter );
      printf0("solve time: %le seconds\n", t1-t0 );
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += iter/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_ITER] += (((double)g.coarse_iter_count)/((double)iter))/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_TIME] += g.coarse_time/((double)g.vt.average_over);
    END_LOCKED_MASTER(threading)
    }
  }
  if ( l->depth == 0 && ( p->timing || p->print ) && !(g.vt.p_end != NULL )  ) {
    START_MASTER(threading)
    if ( g.method != 6 ) prof_print( l );
    END_MASTER(threading)
  }
  
  endProfilingRange(profilingRangeResFgmres);
  if ( g.my_rank == 0 && p->print_iters == 1 ) printf("SOLVER ITERS (depth=%d) : %d\n", l->depth, iter);
  return iter;
}


void bicgstab_PRECISION( gmres_PRECISION_struct *ps, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Uses BiCGstab to solve the system D x = b, where b is taken from ps->b and x is 
* stored in ps->x.                                                              
*********************************************************************************/
  
  vector_PRECISION x, b, r, r_tilde, p, pp, v, s, t; // Krylov subspace size: 5
  complex_PRECISION alpha=1, beta=1, rho=1, rho_old=1, omega=1;
  int iter=0, maxiter;
  double tol, b_norm, r_norm, s_norm;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(ps->v_start, ps->v_end, &start, &end, l, threading);
  
  tol = (l->level==0 && g.num_levels > 1 && g.interpolation )?g.coarse_tol:g.bicgstab_tol;
  maxiter = 1000000; r = ps->r; b = ps->b; x = ps->x; p = ps->w;
  pp = ps->V[0]; r_tilde = ps->V[1]; v = ps->V[2]; s = ps->V[3]; t = ps->V[4];
  
  vector_PRECISION_copy( r, b, start, end, l );
  vector_PRECISION_copy( r_tilde, b, start, end, l );
  vector_PRECISION_define( x, 0, start, end, l );
  vector_PRECISION_define( v, 0, start, end, l );
  vector_PRECISION_define( s, 0, start, end, l );
  vector_PRECISION_define( t, 0, start, end, l );
  b_norm = global_norm_PRECISION( b, ps->v_start, ps->v_end, l, threading );
  r_norm = b_norm;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  START_MASTER(threading)
  printf0("+----------------------------------------------------------+\n");
  END_MASTER(threading)
#endif
  while ( r_norm/b_norm > tol && iter < maxiter ) {
    iter++;
    
    rho_old = rho;
    rho = global_inner_product_PRECISION( r_tilde, r, ps->v_start, ps->v_end, l, threading );
    if ( rho == 0 ) {
      START_MASTER(threading)
      printf0("rho = 0: BiCGstab did not converge.\n");
      END_MASTER(threading)
      break;
    }
    
    if ( iter == 1 ) {
      vector_PRECISION_copy( p, r, start, end, l );
    } else {
      beta = (rho/rho_old)*(alpha/omega);
      vector_PRECISION_saxpy( pp, p,  v, -omega, start, end, l );
      vector_PRECISION_saxpy( p,  r, pp,   beta, start, end, l );
    }    
    apply_operator_PRECISION( v, p, ps, l, threading );
    alpha = rho / global_inner_product_PRECISION( r_tilde, v, ps->v_start, ps->v_end, l, threading );
    vector_PRECISION_saxpy( s, r, v, -alpha, start, end, l );
    s_norm = global_norm_PRECISION( s, ps->v_start, ps->v_end, l, threading );
    
    if ( s_norm/b_norm < tol ) {
      vector_PRECISION_saxpy( x, x, p, alpha, start, end, l );
      break;
    }
    
    apply_operator_PRECISION( t, s, ps, l, threading );
    omega = global_inner_product_PRECISION( t, s, ps->v_start, ps->v_end, l, threading )
          / global_inner_product_PRECISION( t, t, ps->v_start, ps->v_end, l, threading );
    vector_PRECISION_saxpy( x, x, p,  alpha, start, end, l );
    vector_PRECISION_saxpy( x, x, s,  omega, start, end, l );
    vector_PRECISION_saxpy( r, s, t, -omega, start, end, l );
    
    r_norm = global_norm_PRECISION( r, ps->v_start, ps->v_end, l, threading );
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    START_MASTER(threading)
    if ( iter % 100 == 0 ) printf0("| biCGstab relres: %12.6le,  iterations: %-8d     |\n", r_norm/b_norm, iter );
    END_MASTER(threading)
#endif
  }
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  START_MASTER(threading)
  printf0("| biCGstab relres: %12.6le,  iterations: %-8d     |\n", r_norm/b_norm, iter );
  printf0("+----------------------------------------------------------+\n");
  END_MASTER(threading)
#endif
}


void cgn_PRECISION( gmres_PRECISION_struct *ps, level_struct *l, struct Thread *threading ) {
  
/*********************************************************************************
* Uses CGN to solve the system D x = b, where b is taken from ps->b and x is 
* stored in ps->x.                                                              
*********************************************************************************/

  vector_PRECISION r_old, r_new, r_true, p, pp, Dp, x, b;
  complex_PRECISION alpha, beta=0, gamma;
  int maxiter, iter=0;
  double tol, r0_norm, r_norm, prod_rr_old, t0=0, t1=0;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  
  b = ps->b; x = ps->x;
  r_old = ps->V[2]; r_new = ps->V[3]; r_true = ps->r;
  p = ps->w; pp = ps->V[0]; Dp = ps->V[1];
  tol = (l->level==0 && g.num_levels > 1 && g.interpolation )?g.coarse_tol:ps->tol;
  maxiter = ps->num_restart;
  
  START_MASTER(threading)
  if ( ps->timing || ps->print ) t0 = MPI_Wtime();
  END_MASTER(threading)

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(ps->v_start, ps->v_end, &start, &end, l, threading);

  vector_PRECISION_define( x, 0, start, end, l );
  apply_operator_PRECISION( Dp, x, ps, l, threading );
  vector_PRECISION_minus( pp, b, Dp, start, end, l );
  apply_operator_dagger_PRECISION( r_old, pp, ps, l, threading );
  
  vector_PRECISION_copy( p, r_old, start, end, l );
  r0_norm = creal(global_norm_PRECISION( r_old, ps->v_start, ps->v_end, l, threading ));
  prod_rr_old = global_inner_product_PRECISION( r_old, r_old, ps->v_start, ps->v_end, l, threading );
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("\n+----------------------------------------------------------+\n");
    END_MASTER(threading)
  }  
#endif
  while ( sqrt(prod_rr_old) / r0_norm > tol && iter < maxiter ) {
    iter++;
    
    apply_operator_PRECISION( pp, p, ps, l, threading );
    apply_operator_dagger_PRECISION( Dp, pp, ps, l, threading );
    
    gamma = global_inner_product_PRECISION( p, Dp, ps->v_start, ps->v_end, l, threading );
    alpha = prod_rr_old / gamma;
    vector_PRECISION_saxpy( x, x, p, alpha, start, end, l );
    vector_PRECISION_saxpy( r_new, r_old, Dp, -alpha, start, end, l );
    
    gamma = global_inner_product_PRECISION( r_new, r_new, ps->v_start, ps->v_end, l, threading );
    beta = gamma / prod_rr_old;
    
    vector_PRECISION_saxpy( p, r_new, p, beta, start, end, l );
    vector_PRECISION_copy( r_old, r_new, start, end, l );
    prod_rr_old = gamma;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)    
    if ( iter%100 == 0 && ps->print >=1 ) {
      START_MASTER(threading)
      printf0("|      NE rel. res. after  %-6d iterations: %e |\n", iter, sqrt(prod_rr_old)/r0_norm );
      END_MASTER(threading)
    }
#endif
  }
  
  r0_norm = creal(global_norm_PRECISION( b, ps->v_start, ps->v_end, l, threading ));
  apply_operator_PRECISION( Dp, x, ps, l, threading );
  vector_PRECISION_minus( r_true, b, Dp, start, end, l );
  r_norm = creal(global_norm_PRECISION( r_true, ps->v_start, ps->v_end, l, threading ));
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("+----------------------------------------------------------+\n");
    printf0("| switching to CGNR, iter  %-6d true r res: %e |\n", iter, r_norm/r0_norm );
    printf0("+----------------------------------------------------------+\n");
    END_MASTER(threading)
  }
#endif
  
  while ( r_norm / r0_norm > tol && iter < maxiter ) {
    iter++;
    
    apply_operator_PRECISION( pp, p, ps, l, threading );
    apply_operator_dagger_PRECISION( Dp, pp, ps, l, threading );
    
    gamma = global_inner_product_PRECISION( p, Dp, ps->v_start, ps->v_end, l, threading );
    alpha = prod_rr_old / gamma;
    vector_PRECISION_saxpy( x, x, p, alpha, start, end, l );
    vector_PRECISION_saxpy( r_new, r_old, Dp, -alpha, start, end, l );
    
    // residual update
    vector_PRECISION_saxpy( r_true, r_true, pp, -alpha, start, end, l );
    r_norm = creal(global_norm_PRECISION( r_true, ps->v_start, ps->v_end, l, threading ));
    gamma = global_inner_product_PRECISION( r_new, r_new, ps->v_start, ps->v_end, l, threading );
    beta = gamma / prod_rr_old;
    
    vector_PRECISION_saxpy( p, r_new, p, beta, start, end, l );
    vector_PRECISION_copy( r_old, r_new, start, end, l );
    prod_rr_old = gamma;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)    
    if ( iter%100 ==  0 && ps->print >=1 ) {
      START_MASTER(threading)
      printf0("|         rel. res. after  %-6d iterations: %e |\n", iter, r_norm/r0_norm );
      END_MASTER(threading)
    }
#endif
  }
  
  if ( ps->timing || ps->print ) t1 = MPI_Wtime();
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("+----------------------------------------------------------+\n");
    printf0("|          CGN iterations: %-6d                          |\n", iter );
    END_MASTER(threading)
    apply_operator_PRECISION( Dp, x, ps, l, threading );
    vector_PRECISION_minus( pp, b, Dp, start, end, l );
    beta = creal(global_norm_PRECISION( pp, ps->v_start, ps->v_end, l, threading ));
    START_MASTER(threading)
    if ( ps->timing ) printf0("| exact relative residual: ||r||/||b|| = %e      |\n", creal(beta/r0_norm) );
    printf0("| elapsed wall clock time: %-12g seconds            |\n", t1-t0 );
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
#ifdef CUDA_OPT
    printf0("|        max used mem/GPU: %-8.2le GB                     |\n", g.max_gpu_storage/1024.0 );
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }
  
  START_LOCKED_MASTER(threading)
  if ( l->level == 0 )
    g.coarse_iter_count += iter;

  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += (iter)/((double)g.vt.average_over);
    }
  }
  END_LOCKED_MASTER(threading)
}


int arnoldi_step_PRECISION( vector_PRECISION *V, vector_PRECISION *Z, vector_PRECISION w,
                            complex_PRECISION **H, complex_PRECISION* buffer, int j, void (*prec)(),
                            complex_PRECISION shift, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {
  RangeHandleType profilingRangeStep = startProfilingRange("Arnoldi step (PRECISION)");

/*********************************************************************************
* Extends the Arnoldi basis by one vector.
* - vector_PRECISION *V: Contains the Arnoldi basis vectors.
* - vector_PRECISION *Z: If a right precond. P is used, contains P*V[j] for all j.
* - vector_PRECISION w: Will be appended to existing Arnoldi basis at 
*   position j+1.
* - complex_PRECISION **H: Contains full Hessenberg matrix from the Arnoldi 
*   decomposition (columnmajor!)
* - complex_PRECISION* buffer: Buffer for local inner products.
* - int j: index of the new Arnoldi vector to be orthonormalized
*   against all previous ones.
* - void (*prec)(): Function pointer to preconditioner (can be NULL if no 
*   preconditioning is used).
* - complex_PRECISION shift: Denotes the dirac shift (can be 0).
*********************************************************************************/
//#ifdef SINGLE_ALLREDUCE_ARNOLDI
#if 0
#ifdef PIPELINED_ARNOLDI
  if ( l->level == 0 && l->depth > 0 ) {
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    MPI_Request req;
    MPI_Status stat;
    int start, end, i;
    const complex_PRECISION sigma = 0;
    compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
    
    if ( j == 0 )
      vector_PRECISION_copy( Z[0], V[0], start, end, l );
    else
      vector_PRECISION_copy( V[j], Z[j], start, end, l );
  
    complex_PRECISION tmp[j+1];
    process_multi_inner_product_PRECISION( j+1, tmp, V, V[j], p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    for( i=0; i<=j; i++ ) {
      buffer[i] = tmp[i];
    }
    if ( g.num_processes > 1 ) {
      MPI_Iallreduce( buffer, H[MAX(0,j-1)], j+1, MPI_COMPLEX_PRECISION, MPI_SUM,
                      (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm, &req );
    } else {
      for( i=0; i<=j; i++ )
        H[MAX(0,j-1)][i] = buffer[i];
    }
    PROF_PRECISION_STOP( _ALLR, 1 );
    END_MASTER(threading)
    apply_operator_PRECISION( Z[j+1], Z[j], p, l, threading );

    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    if ( g.num_processes > 1 ) {
      MPI_Wait( &req, &stat );
    }
    PROF_PRECISION_STOP( _ALLR, 0 );
    if ( j > 0 ) {
      for ( i=0; i<j; i++ )
        H[j-1][j] -= conj( H[j-1][i] )*H[j-1][i];
    }
    H[MAX(0,j-1)][j] = sqrt( creal( H[MAX(0,j-1)][j] ) );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading) 
    
    for( i=0; i<j; i++ )
      vector_PRECISION_saxpy( V[j], V[j], V[i], -H[j-1][i], start, end, l );
    
    vector_PRECISION_real_scale( V[j], V[j], 1/H[MAX(0,j-1)][j], start, end, l );
    
    START_MASTER(threading)
    if ( j > 0 ) {
      H[j-1][j-1] += sigma;
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if ( j == 0 ) {
      if ( sigma ) vector_PRECISION_saxpy( Z[j+1], Z[j+1], Z[j], -sigma, start, end, l );
    } else {
      for( i=0; i<j; i++ )
        vector_PRECISION_saxpy( Z[j+1], Z[j+1], Z[i+1], -H[j-1][i], start, end, l );
    }
    
    vector_PRECISION_real_scale( Z[j+1], Z[j+1], 1/H[MAX(0,j-1)][j], start, end, l );
    
  } else {
#endif
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    int start, end, i;
    const complex_PRECISION sigma = 0;
    compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
    
    if ( prec != NULL ) {
      if ( p->kind == _LEFT ) {
        apply_operator_PRECISION( Z[0], V[j], p, l, threading );
        if ( shift ) vector_PRECISION_saxpy( Z[0], Z[0], V[j], shift, start, end, l );
        prec( V[j+1], NULL, Z[0], _NO_RES, l, threading );
        if ( sigma ) vector_PRECISION_saxpy( V[j+1], V[j+1], V[j], -sigma, start, end, l );
      } else {
        if ( l->level == 0 ) {
          prec( Z[j], NULL, V[j], _NO_RES, l, threading );
          apply_operator_PRECISION( V[j+1], Z[j], p, l, threading );
        } else {
          if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2 ) ) {
            prec( Z[j], V[j+1], V[j], _NO_RES, l, threading );
            // obtains w = D * Z[j] from Schwarz
          } else {
            prec( Z[j], NULL, V[j], _NO_RES, l, threading );
            apply_operator_PRECISION( V[j+1], Z[j], p, l, threading ); // w = D*Z[j]
          }
        }
        if ( shift ) vector_PRECISION_saxpy( V[j+1], V[j+1], Z[j], shift, start, end, l );
        if ( sigma ) vector_PRECISION_saxpy( V[j+1], V[j+1], V[j], -sigma, start, end, l );
      }
    } else {
      apply_operator_PRECISION( V[j+1], V[j], p, l, threading ); // w = D*V[j]
      if ( shift-sigma ) vector_PRECISION_saxpy( V[j+1], V[j+1], V[j], shift-sigma, start, end, l );
    }
    
    complex_PRECISION tmp[j+2];
    process_multi_inner_product_PRECISION( j+2, tmp, V, V[j+1], p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    for( i=0; i<=j+1; i++ ) {
      buffer[i] = tmp[i];
    }
    
    if ( g.num_processes > 1 ) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, H[j], j+2, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    } else {
      for( i=0; i<=j+1; i++ )
        H[j][i] = buffer[i];
    }  
    for ( i=0; i<=j; i++ )
      H[j][j+1] -= conj( H[j][i] )*H[j][i];
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    if ( creal( H[j][j+1] ) < 0 )
      return 0;
    START_MASTER(threading)
    H[j][j+1] = sqrt( creal( H[j][j+1] ) );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    for( i=0; i<=j; i++ )
      vector_PRECISION_saxpy( V[j+1], V[j+1], V[i], -H[j][i], start, end, l );
    
    vector_PRECISION_real_scale( V[j+1], V[j+1], 1/H[j][j+1], start, end, l );
    START_LOCKED_MASTER(threading)
    H[j][j] += sigma;
    END_LOCKED_MASTER(threading)
#ifdef PIPELINED_ARNOLDI
  }
#endif
#else
  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
  int i;
  // start and end indices for vector functions depending on thread
  int start, end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  
  if ( prec != NULL ) {
    if ( p->kind == _LEFT ) {
      apply_operator_PRECISION( Z[0], V[j], p, l, threading );
      if ( shift ) vector_PRECISION_saxpy( Z[0], Z[0], V[j], shift, start, end, l );
      prec( w, NULL, Z[0], _NO_RES, l, threading );
    } else {
      if ( l->level == 0 ) { 
        prec( Z[j], NULL, V[j], _NO_RES, l, threading );
        apply_operator_PRECISION( w, Z[j], p, l, threading );
        //apply_operator_PRECISION( w, Z[j], p, l, threading );
      } else {
        if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2 ) ) {
          prec( Z[j], w, V[j], _NO_RES, l, threading );
          // obtains w = D * Z[j] from Schwarz
        } else {
          prec( Z[j], NULL, V[j], _NO_RES, l, threading );
          apply_operator_PRECISION( w, Z[j], p, l, threading ); // w = D*Z[j]
        }
        if ( shift ) vector_PRECISION_saxpy( w, w, Z[j], shift, start, end, l );
      }
    }
  } else {
    apply_operator_PRECISION( w, V[j], p, l, threading ); // w = D*V[j
    if ( shift ) vector_PRECISION_saxpy( w, w, V[j], shift, start, end, l );
  }

// LAST STAGE
#ifdef GCRODR
//#if 0
  // orthogonalize against Ck whenever necessary
  if ( l->level==0 && p->gcrodr_PRECISION.orth_against_Ck == 1 ) {
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    int k = p->gcrodr_PRECISION.k;
    vector_PRECISION *Ck = p->gcrodr_PRECISION.C;
    complex_PRECISION **B = p->gcrodr_PRECISION.ort_B;
    // buffer
    complex_PRECISION *bf = p->gcrodr_PRECISION.Bbuff[0];

    complex_PRECISION tmpx[k+j+2];

    // merging Ck and V into a single vector of pointers
    complex_PRECISION* VCk[k+j+2];
    for( i=0;i<k;i++ ){ VCk[i] = Ck[i]; }
    for( i=0;i<(j+1);i++ ){ VCk[k+i] = V[i]; }
    VCk[k+j+1] = w; 

    process_multi_inner_product_PRECISION( k+j+2, tmpx, VCk, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    // buffer is of length m, and k<m
    for ( i=0; i<(k+j+2); i++ )
      buffer[i] = tmpx[i];
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, bf, k+j+2, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    } else {
      for( i=0; i<(k+j+2); i++ )
        bf[i] = buffer[i];
    }

    // copy the B coefficients to the corresponding matrix
    memcpy( B[j], bf, sizeof(complex_PRECISION)*k );
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    for( i=0; i<k; i++ )
      vector_PRECISION_saxpy( w, w, Ck[i], -B[j][i], start, end, l );

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_MASTER(threading)
    // copy the H coefficients to the corresponding matrix
    memcpy( H[j], bf+k, sizeof(complex_PRECISION)*(j+2) );
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    for( i=0; i<=j; i++ )
      vector_PRECISION_saxpy( w, w, V[i], -H[j][i], start, end, l );
#ifdef REORTH
    // re-orthogonalization
    process_multi_inner_product_PRECISION( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    for( i=0; i<=j; i++ )
      buffer[i] = tmp[i];
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, tmp, j+1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    }

    for( i=0; i<=j; i++ )
      H[j][i] += tmp[i];

    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    for( i=0; i<=j; i++ )
      vector_PRECISION_saxpy( w, w, V[i], -tmp[i], start, end, l );
#endif

  } else {

    // orthogonalization
    complex_PRECISION tmp[j+2];
    complex_PRECISION *V_buff[j+2];

    if ( l->level==0 ) {
      for (i=0; i < j+1; i++) V_buff[i] = V[i];
      V_buff[j+1] = w;
    }

    if ( l->level==0 ) process_multi_inner_product_PRECISION( j+2, tmp, V_buff, w, p->v_start, p->v_end, l, threading );
    else process_multi_inner_product_PRECISION( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    if ( l->level==0 ) {
      for( i=0; i<=j+1; i++ )
        buffer[i] = tmp[i];
    } else {
      for( i=0; i<=j; i++ )
        buffer[i] = tmp[i];
    }
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_PRECISION_START( _ALLR );
      if ( l->level==0 ) MPI_Allreduce( buffer, H[j], j+2, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      else MPI_Allreduce( buffer, H[j], j+1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    } else {
      if ( l->level==0 ) {
        for( i=0; i<=j+1; i++ )
          H[j][i] = buffer[i];
      } else {
        for( i=0; i<=j; i++ )
          H[j][i] = buffer[i];
      }
    }
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    for( i=0; i<=j; i++ )
      vector_PRECISION_saxpy( w, w, V[i], -H[j][i], start, end, l );
#ifdef REORTH
    // re-orthogonalization
    process_multi_inner_product_PRECISION( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    for( i=0; i<=j; i++ )
      buffer[i] = tmp[i];
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, tmp, j+1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    }

    for( i=0; i<=j; i++ )
      H[j][i] += tmp[i];

    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    for( i=0; i<=j; i++ )
      vector_PRECISION_saxpy( w, w, V[i], -tmp[i], start, end, l );
#endif

  }

#else

  // orthogonalization
  complex_PRECISION tmp[j+2];
  complex_PRECISION *V_buff[j+2];

  if ( l->level==0 ) {
    for (i=0; i < j+1; i++) V_buff[i] = V[i];
    V_buff[j+1] = w;
  }

  if ( l->level==0 ) process_multi_inner_product_PRECISION( j+2, tmp, V_buff, w, p->v_start, p->v_end, l, threading );
  else process_multi_inner_product_PRECISION( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  if ( l->level==0 ) {
    for( i=0; i<=j+1; i++ )
      buffer[i] = tmp[i];
  } else {
    for( i=0; i<=j; i++ )
      buffer[i] = tmp[i];
  }
  if ( g.num_processes > 1 ) {
    //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
    PROF_PRECISION_START( _ALLR );
    if ( l->level==0 ) MPI_Allreduce( buffer, H[j], j+2, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    else MPI_Allreduce( buffer, H[j], j+1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
  } else {
    if ( l->level==0 ) {
      for( i=0; i<=j+1; i++ )
        H[j][i] = buffer[i];
    } else {
      for( i=0; i<=j; i++ )
        H[j][i] = buffer[i];
    }
  }
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  for( i=0; i<=j; i++ )
    vector_PRECISION_saxpy( w, w, V[i], -H[j][i], start, end, l );
#ifdef REORTH
  // re-orthogonalization
  process_multi_inner_product_PRECISION( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  for( i=0; i<=j; i++ )
    buffer[i] = tmp[i];
  if ( g.num_processes > 1 ) {
    //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( buffer, tmp, j+1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
  }
  
  for( i=0; i<=j; i++ )
    H[j][i] += tmp[i];

  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  for( i=0; i<=j; i++ )
    vector_PRECISION_saxpy( w, w, V[i], -tmp[i], start, end, l );
#endif

#endif // from GCRO-DR
  
  //// normalization
  //complex_PRECISION tmp2 = global_norm_PRECISION( w, p->v_start, p->v_end, l, threading );
  //START_MASTER(threading)
  //H[j][j+1] = tmp2;
  //END_MASTER(threading)
  //SYNC_MASTER_TO_ALL(threading)

//=============================================================
  if ( l->level==0 ) {
    START_MASTER(threading)
    complex_PRECISION tmp = H[j][j+1];
    for ( i=0; i<=j; i++ )
      tmp -= conj_PRECISION( H[j][i] )*H[j][i];

// LAST STAGE
#ifdef GCRODR
//#if 0
    if ( l->level==0 && p->gcrodr_PRECISION.orth_against_Ck == 1 )
    {      
      int k = p->gcrodr_PRECISION.k;
      complex_PRECISION **B = p->gcrodr_PRECISION.ort_B;

      for( i=0; i<k; i++ )
        tmp -= conj_PRECISION( B[j][i] )*B[j][i];
    }
 #endif
      
    H[j][j+1] = tmp;
    ((complex_PRECISION*)threading->workspace)[0] = creal_PRECISION(tmp);
    END_MASTER(threading)
      
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    if ( creal_PRECISION(((complex_PRECISION*)threading->workspace)[0]) < 0.)
    {
      PRECISION tmp3 = global_norm_PRECISION( w, p->v_start, p->v_end, l, threading );
      START_MASTER(threading)
      H[j][j+1] = tmp3;
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    else
    {
      START_MASTER(threading)
      H[j][j+1] = sqrt(creal_PRECISION(H[j][j+1]));
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
    }

  } else {
    PRECISION tmp2 = global_norm_PRECISION( w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    H[j][j+1] = tmp2;
    END_MASTER(threading)
  }
//=============================================================

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)
  
  // V_j+1 = w / H_j+1,j
  if ( cabs_PRECISION( H[j][j+1] ) > 1e-15 )
    vector_PRECISION_real_scale( V[j+1], w, 1/H[j][j+1], start, end, l );

  START_MASTER(threading)

#if defined(GCRODR) || defined(POLYPREC)
//#if defined(POLYPREC)
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)

  //int jx = j-1;
  int jx;
  if ( j==0 ) {
    jx = j;
  } else {
    jx = j-1;
  }
#else
  int jx = j;
#endif
#endif

  // copy of Hesselnberg matrix (only level=0 currently)
#if defined(GCRODR) && defined(POLYPREC)
//#if 0
  if (l->dup_H==1 && l->level==0)
  {
    memcpy( p->gcrodr_PRECISION.eigslvr.Hc[jx], H[jx], sizeof(complex_PRECISION)*(jx+2) );
    memset( p->gcrodr_PRECISION.eigslvr.Hc[jx]+jx+2, 0.0, sizeof(complex_PRECISION)*(p->restart_length + 1 - (jx+2)) );
  }
#elif defined(GCRODR)
//#elif 0
  if (l->dup_H==1 && l->level==0)
  {
    memcpy( p->gcrodr_PRECISION.eigslvr.Hc[jx], H[jx], sizeof(complex_PRECISION)*(jx+2) );
    memset( p->gcrodr_PRECISION.eigslvr.Hc[jx]+jx+2, 0.0, sizeof(complex_PRECISION)*(p->restart_length + 1 - (jx+2)) );
  }
#elif defined(POLYPREC)
  if (l->dup_H==1 && l->level==0)
  {
    memcpy( p->polyprec_PRECISION.eigslvr.Hc[jx], H[jx], sizeof(complex_PRECISION)*(jx+2) );
    memset( p->polyprec_PRECISION.eigslvr.Hc[jx]+jx+2, 0.0, sizeof(complex_PRECISION)*(p->restart_length + 1 - (jx+2)) );
  }
#endif

  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

#endif
  endProfilingRange(profilingRangeStep);
  return 1;
}


void qr_update_PRECISION( complex_PRECISION **H, complex_PRECISION *s,
                          complex_PRECISION *c, complex_PRECISION *gamma, int j,
                          level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Applies one Givens rotation to the Hessenberg matrix H in order to solve the 
* least squares problem in (F)GMRES for computing the solution.
* - complex_PRECISION **H: Hessenberg matrix from Arnoldi decomposition
* - complex_PRECISION *s: sin values from givens rotations
* - complex_PRECISION *c: cos valies from givens rotations
* - complex_PRECISION *gamma: Approximation to residual from every step
* - int j: Denotes current iteration.
*********************************************************************************/  

  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  START_MASTER(threading)
  
  PROF_PRECISION_START( _SMALL1 );
  
  int i;
  complex_PRECISION beta;
  
  // update QR factorization
  // apply previous Givens rotation
  for ( i=0; i<j; i++ ) {
    beta = (-s[i])*H[j][i] + (c[i])*H[j][i+1];
    H[j][i] = conj_PRECISION(c[i])*H[j][i] + conj_PRECISION(s[i])*H[j][i+1];
    H[j][i+1] = beta;
  }
  // compute current Givens rotation
  beta = (complex_PRECISION) sqrt( NORM_SQUARE_PRECISION(H[j][j]) + NORM_SQUARE_PRECISION(H[j][j+1]) );
  s[j] = H[j][j+1]/beta; c[j] = H[j][j]/beta;
  // update right column
  gamma[j+1] = (-s[j])*gamma[j]; gamma[j] = conj_PRECISION(c[j])*gamma[j];
  // apply current Givens rotation
  H[j][j] = beta; H[j][j+1] = 0;
  
  PROF_PRECISION_STOP( _SMALL1, 6*j+6 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
}


void compute_solution_PRECISION( vector_PRECISION x, vector_PRECISION *V, complex_PRECISION *y,
                                 complex_PRECISION *gamma, complex_PRECISION **H, int j, int ol,
                                 gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {
  RangeHandleType profilingRangeSolution = startProfilingRange("compute solution (PRECISION)");
  
  int i, k;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  START_MASTER(threading)
  
  PROF_PRECISION_START( _SMALL2 );
  
  // backward substitution
  for ( i=j; i>=0; i-- ) {
    y[i] = gamma[i];
    for ( k=i+1; k<=j; k++ ) {
      y[i] -= H[k][i]*y[k];
    }
    y[i] /= H[i][i];
  }
  
  PROF_PRECISION_STOP( _SMALL2, ((j+1)*(j+2))/2 + j+1 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // x = x + V*y
  if ( ol ) {
    for ( i=0; i<=j; i++ )
      vector_PRECISION_saxpy( x, x, V[i], y[i], start, end, l );
  } else {
    vector_PRECISION_scale( x, V[0], y[0], start, end, l );
    for ( i=1; i<=j; i++ )
      vector_PRECISION_saxpy( x, x, V[i], y[i], start, end, l );
  }
  endProfilingRange(profilingRangeSolution);
}


void local_minres_PRECISION( vector_PRECISION phi, vector_PRECISION eta, vector_PRECISION latest_iter,
                             int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
/*********************************************************************************
* Minimal Residual iteration solver used to solve the block systems
*     blockD phi = eta
* within the Schwarz method, phi contains an initial guess and its updated version
* is returned after the block solve has been performed.
* eta is overwritten by the block residual r.
* To calculate the missing contributions to r on the current Schwarz block
* coming from outside of the block, an update "phi_new - phi_old" is returned in
* latest_iter -> cheaper residual update in the Schwarz method
*********************************************************************************/
  
  START_UNTHREADED_FUNCTION(threading)

  int i, end = (g.odd_even&&l->depth==0)?start+12*s->num_block_even_sites:start+s->block_vector_size,
      n = l->block_iter;
  vector_PRECISION Dr = s->local_minres_buffer[0];
  vector_PRECISION r = s->local_minres_buffer[1];
  vector_PRECISION lphi = s->local_minres_buffer[2];
  complex_PRECISION alpha;
  void (*block_op)() = (l->depth==0)?(g.odd_even?apply_block_schur_complement_PRECISION:block_d_plus_clover_PRECISION)
                                    :coarse_block_operator_PRECISION;
  
  vector_PRECISION_copy( r, eta, start, end, l );
  vector_PRECISION_define( lphi, 0, start, end, l );
  
  for ( i=0; i<n; i++ ) {

    // Dr = blockD*r
    block_op( Dr, r, start, s, l, no_threading );

    // alpha = <Dr,r>/<Dr,Dr>
    alpha = local_xy_over_xx_PRECISION( Dr, r, start, end, l );

    // phi += alpha * r
    vector_PRECISION_saxpy( lphi, lphi, r, alpha, start, end, l );
    // r -= alpha * Dr
    vector_PRECISION_saxpy( r, r, Dr, -alpha, start, end, l );
  }
  
  if ( latest_iter != NULL ) vector_PRECISION_copy( latest_iter, lphi, start, end, l );
  if ( phi != NULL ) vector_PRECISION_plus( phi, phi, lphi, start, end, l );
  vector_PRECISION_copy( eta, r, start, end, l );

  END_UNTHREADED_FUNCTION(threading)
}


// used as smoother at the moment
int fgcr_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) { 

/*********************************************************************************
* Uses FGCR to solve the system D x = b, where b is taken from p->b and x is 
* stored in p->x.                                                              
*********************************************************************************/

  int k, j, iter=0, ol, start, end;

  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );

  START_MASTER(threading)
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  if ( p->print ) printf0("+----------------------------------------------------------+\n");
#endif
  END_MASTER(threading)

  if ( p->initial_guess_zero==1 ) {
    vector_PRECISION_define( p->x, 0, start, end, l );
  }

  // in the notation of the ORISE paper, p->V are the p vectors,
  // and p->Z the z vectors. Let us make it simpler then by creating
  // two new pointers
  vector_PRECISION *Z = p->Z;
  vector_PRECISION *P = p->V;

  complex_PRECISION *buffer = p->gcr_buffer_dotprods;
  complex_PRECISION *betas  = p->gcr_betas_dotprods;

  for ( ol=0;ol<p->num_restart;ol++ ) {

    // compute the residual
    if ( p->initial_guess_zero==1 && ol==0 ) {
      // if the initial guess is zero and we're at the first iteration,
      // just equal the residual to the RHS
      vector_PRECISION_copy( p->r, p->b, start, end, l );
    } else {
      apply_operator_PRECISION( p->w, p->x, p, l, threading );
      vector_PRECISION_minus( p->r, p->b, p->w, start, end, l );
    }

    // copy r into P[0]
    vector_PRECISION_copy( P[0], p->r, start, end, l );
    // and build the first Z vector
    apply_operator_PRECISION( Z[0], P[0], p, l, threading );

    complex_PRECISION alpha;
    complex_PRECISION deltas[p->restart_length];

    for ( k=0;k<p->restart_length;k++ ) {
      // global number of iters
      iter++;

      deltas[k] = global_inner_product_PRECISION( Z[k], Z[k], p->v_start, p->v_end, l, threading );
      alpha     = global_inner_product_PRECISION( Z[k], p->r, p->v_start, p->v_end, l, threading ) / deltas[k];

      vector_PRECISION_saxpy( p->x, p->x, P[k], alpha,  start, end, l );
      vector_PRECISION_saxpy( p->r, p->r, Z[k], -alpha, start, end, l );

      apply_operator_PRECISION( p->w, p->r, p, l, threading );

      {
        // orthogonalization
        complex_PRECISION tmp[k+1];

        process_multi_inner_product_PRECISION( k+1, tmp, Z, p->w, p->v_start, p->v_end, l, threading );

        START_MASTER(threading)
        for ( j=0; j<=k; j++ )
          buffer[j] = tmp[j];

        if ( g.num_processes > 1 ) {
          PROF_PRECISION_START( _ALLR );
          MPI_Allreduce( buffer, betas, k+1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
          PROF_PRECISION_STOP( _ALLR, 1 );
        } else {
          for( j=0; j<=k; j++ )
            betas[j] = buffer[j];
        }

        for ( j=0; j<=k; j++ )
          betas[j] = -betas[j]/deltas[j];
        END_MASTER(threading)
        SYNC_CORES(threading)
      }

      //for ( j=0;j<=k;j++ ) {
      //  betas[j] = -global_inner_product_PRECISION( Z[j], p->w, p->v_start, p->v_end, l, threading ) / deltas[j];
      //}

      vector_PRECISION_copy( P[k+1], p->r, start, end, l );
      for ( j=0;j<=k;j++ ) { vector_PRECISION_saxpy( P[k+1], P[k+1], P[j], betas[j],  start, end, l ); }
      vector_PRECISION_copy( Z[k+1], p->w, start, end, l );
      for ( j=0;j<=k;j++ ) { vector_PRECISION_saxpy( Z[k+1], Z[k+1], Z[j], betas[j],  start, end, l ); }
    }
  }

  return iter;
}

#ifdef RICHARDSON_SMOOTHER
void richardson_update_omega_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  // BPI = block power iteration

  int start, end, i, j, k;
  PRECISION norm;
  complex_PRECISION lmaxb, dot_prod;

  // if e.g. we want 2 shifts in Richardson, we use 2 BPI vectors .. but this code
  // is in principle ready to apply block power iteration in a more general way
  int bpi_base = 1 * p->richardson_sub_degree;

  // allocate the BPI base
  vector_PRECISION *V = NULL;
  PUBLIC_MALLOC( V, vector_PRECISION, bpi_base );
  V[0] = NULL;
  PUBLIC_MALLOC( V[0], complex_PRECISION, l->inner_vector_size * bpi_base );
  for( i=1;i<bpi_base;i++ ){
    V[i] = V[0] + i*l->inner_vector_size;
  }

  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );

  // number of power iteration iters, just rough computation
  int bpi_iters = g.smoother_richardson_BPI_iters;

  // do a bit of power iteration to roughly estimate the max eigenvalue

  // initially set the BPI vectors to random
  START_MASTER(threading)
  for ( i=0;i<bpi_base;i++ ) {
    vector_PRECISION_define_random( V[i], 0, l->inner_vector_size, l );
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // and then normalize the BPI vectors
  for ( i=0;i<bpi_base;i++ ) {
    norm = global_norm_PRECISION( V[i], p->v_start, p->v_end, l, threading );
    vector_PRECISION_scale( V[i], V[i], 1.0/norm, start, end, l );
  }

  // do BPI

  for ( i=0;i<bpi_iters;i++ ) {
    // 1. apply D on all the BPI vectors
    for ( j=0;j<bpi_base;j++ ) {
      apply_operator_PRECISION( p->w, V[j], p, l, threading );
      vector_PRECISION_copy( V[j], p->w, start, end, l );
    }
    // 2. orthonormalize the vectors in the base
    for ( j=0;j<bpi_base;j++ ) {
      // first project
      for ( k=0;k<j;k++ ) {
        // dot product
        dot_prod = global_inner_product_PRECISION( V[k], V[j], p->v_start, p->v_end, l, threading );
        // axpy
        vector_PRECISION_saxpy( V[j], V[j], V[k], -dot_prod, start, end, l );
      }
      // then normalize the projected vector
      norm = global_norm_PRECISION( V[j], p->v_start, p->v_end, l, threading );
      vector_PRECISION_scale( V[j], V[j], 1.0/norm, start, end, l );
    }
  }

  // extract the largest Rayleigh quotients

  complex_PRECISION lmax[p->richardson_sub_degree];
  for ( i=0;i<p->richardson_sub_degree;i++ ) {
    lmax[i] = 0.0;
  }
  for ( i=0;i<bpi_base;i++ ) {
    // extract the i-th Rayleigh quotient
    apply_operator_PRECISION( p->w, V[i], p, l, threading );
    norm = global_norm_PRECISION( V[i], p->v_start, p->v_end, l, threading );
    lmaxb = global_inner_product_PRECISION( V[i], p->w, p->v_start, p->v_end, l, threading ) / (norm*norm);
    // if this i-th Rayleigh quotient is larger than the j-th stored one, replace
    for ( j=0;j<p->richardson_sub_degree;j++ ) {
      if ( cabs_PRECISION(lmaxb)>cabs_PRECISION(lmax[j]) ) { lmax[j] = lmaxb; break; }
    }
  }

  //p->omega[i%p->richardson_sub_degree]

  START_MASTER(threading)
  for ( i=0;i<p->richardson_sub_degree;i++ ) {
    p->omega[i] = cabs_PRECISION(lmax[i]);
    //p->omega[i] = 1.0 / (2.0*p->omega[i] / 3.0);
    p->omega[i] = p->richardson_factor / (p->omega[i]);
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  PUBLIC_FREE( V[0], complex_PRECISION, l->inner_vector_size * bpi_base );
  PUBLIC_FREE( V, vector_PRECISION, bpi_base );
}

// used as smoother at the moment
int richardson_PRECISION_cpu( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  if ( p->richardson_update_omega==1 ) {
    richardson_update_omega_PRECISION( p, l, threading );
    START_MASTER(threading)
    p->richardson_update_omega = 0;
    END_MASTER(threading)
  }

  int start, end, i;
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );
  int n = p->num_restart * p->restart_length;

  // initial guess to zero if necessary
  if ( p->initial_guess_zero == _NO_RES ) {
    vector_PRECISION_define( p->x, 0, start, end, l );
  }

  for ( i=0; i<n; i++ ) {
    // 1. compute residual
    if ( i==0 && p->initial_guess_zero==_NO_RES ) {
      vector_PRECISION_copy( p->r, p->b, start, end, l );
    } else {
      apply_operator_PRECISION( p->w, p->x, p, l, threading );
      vector_PRECISION_minus( p->r, p->b, p->w, start, end, l );
    }

    // 2. update solution
    vector_PRECISION_saxpy( p->x, p->x, p->r, p->omega[i%p->richardson_sub_degree], start, end, l );
  }

  return n;
}
#endif
