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
#include "proxies/linsolve_proxy_PRECISION.h"
#include "proxies/operator_proxy_PRECISION.h"

void inv_iter_2lvl_extension_setup_PRECISION( int setup_iter, level_struct *l, struct Thread *threading );
void inv_iter_inv_fcycle_PRECISION( int setup_iter, level_struct *l, struct Thread *threading );
void testvector_analysis_PRECISION( vector_PRECISION *test_vectors, level_struct *l, struct Thread *threading );
void testvector_max_residual_PRECISION( vector_PRECISION *test_vectors, PRECISION *largest_eigen_res, level_struct *l, struct Thread *threading );
void read_tv_from_file_PRECISION( level_struct *l, struct Thread *threading );

void coarse_grid_correction_PRECISION_setup( level_struct *l, struct Thread *threading ) {

  if ( !l->idle ) {

    START_LOCKED_MASTER(threading)
    coarse_operator_PRECISION_alloc( l );
#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
    coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
    END_LOCKED_MASTER(threading)
#else
    END_LOCKED_MASTER(threading)
    coarse_operator_PRECISION_setup_vectorized( l->is_PRECISION.op, l, threading );
#endif

    START_LOCKED_MASTER(threading)
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_PRECISION_alloc( &(l->next_level->s_PRECISION), l->next_level );
        schwarz_layout_PRECISION_define( &(l->next_level->s_PRECISION), l->next_level );
      } else {
        operator_PRECISION_alloc( &(l->next_level->s_PRECISION.op), _ORDINARY, l->next_level );
        operator_PRECISION_define( &(l->next_level->s_PRECISION.op), l->next_level );
        interpolation_PRECISION_alloc( l->next_level );
      }
    } else {
      interpolation_PRECISION_dummy_alloc( l->next_level );
    }

    conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );

    END_LOCKED_MASTER(threading)

    if ( !l->next_level->idle && l->next_level->level > 0 ) {
      START_LOCKED_MASTER(threading)
      schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( g.method >= 4 && g.odd_even ) {
        START_LOCKED_MASTER(threading)
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level );
        END_LOCKED_MASTER(threading)
      }
      coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
      START_LOCKED_MASTER(threading)
      l->next_level->p_PRECISION.op = &(l->next_level->s_PRECISION.op);
      END_LOCKED_MASTER(threading)
    }
    if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
      START_LOCKED_MASTER(threading)
      coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level );
      END_LOCKED_MASTER(threading)
    } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
      coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
    }
  }

  if ( l->next_level->level > 0 ) {
    next_level_setup( NULL, l->next_level, threading );
    START_LOCKED_MASTER(threading)
    if ( !l->next_level->idle )
      interpolation_PRECISION_alloc( l->next_level );
    END_LOCKED_MASTER(threading)
    SYNC_HYPERTHREADS(threading)
    if ( !l->idle ) {
      for ( int i=0; i<MIN(l->next_level->num_eig_vect,l->num_eig_vect); i++ ) {
        restrict_PRECISION( l->next_level->is_PRECISION.test_vector[i], l->is_PRECISION.test_vector[i], l, threading );
      }
      START_LOCKED_MASTER(threading)
      for ( int i=MIN(l->next_level->num_eig_vect,l->num_eig_vect); i<l->next_level->num_eig_vect; i++ ) {
        if ( !l->next_level->idle )
          vector_PRECISION_define_random( l->next_level->is_PRECISION.test_vector[i], 0,
                                          l->next_level->inner_vector_size, l->next_level );
      }
      END_LOCKED_MASTER(threading)
    }
    if ( !l->next_level->idle )
      interpolation_PRECISION_define( NULL, l->next_level, threading );

    coarse_grid_correction_PRECISION_setup( l->next_level, threading );
  }
}


void iterative_PRECISION_setup( int setup_iter, level_struct *l, struct Thread *threading ) {
  if ( l->depth == 0 ) {
    switch ( g.interpolation ) {
      case 2: inv_iter_inv_fcycle_PRECISION( setup_iter, l, threading ); break;
      case 3: inv_iter_inv_fcycle_PRECISION( setup_iter, l, threading ); break;
      case 4: read_tv_from_file_PRECISION( l, threading ); break;
      default: inv_iter_2lvl_extension_setup_PRECISION( setup_iter, l, threading ); break;
    }
  }

  level_struct *lp = l;
  while( lp->level > 0 ) {
    testvector_analysis_PRECISION( lp->is_PRECISION.test_vector, lp, threading );
    lp = lp->next_level;
    if ( lp == NULL )
      break;
  }
}


void read_tv_from_file_PRECISION( level_struct *l, struct Thread *threading ) {

  if ( l->depth == 0 ) {
    if ( g.tv_io_single_file ) {
      START_LOCKED_MASTER(threading)
      vector_io_single_file( NULL, NULL, g.tv_io_file_name, _READ, l->num_eig_vect, "test vectors", l );
      END_LOCKED_MASTER(threading)
      re_setup_PRECISION( l, threading );
    } else {
      START_LOCKED_MASTER(threading)

      int n = l->num_eig_vect, i;
      //char filename[STRINGLENGTH+1];
      char filename[STRINGLENGTH+11];
      vector_double tmp = NULL;

      MALLOC( tmp, complex_double, l->inner_vector_size );

      for ( i=0; i<n; i++ ) {
        sprintf( filename, "%s.%02d", g.tv_io_file_name, i );
        printf0("%s.%02d\n", g.tv_io_file_name, i );
        vector_io( (double*)tmp, filename, _READ, l );
        trans_PRECISION( l->is_PRECISION.test_vector[i], tmp, l->s_PRECISION.op.translation_table, l, no_threading );
      }

      FREE( tmp, complex_double, l->inner_vector_size );

      END_LOCKED_MASTER(threading)

      re_setup_PRECISION( l, threading );
    }
  }
}


void coarse_grid_correction_PRECISION_free( level_struct *l ) {

  next_level_free( l->next_level );

  if ( !l->idle ) {
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_PRECISION_free( &(l->next_level->s_PRECISION), l->next_level );
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_free_PRECISION( l->next_level );
        }
      } else {
        operator_PRECISION_free( &(l->next_level->s_PRECISION.op), _ORDINARY, l->next_level );
        interpolation_PRECISION_free( l->next_level );
        if ( g.odd_even )
          coarse_oddeven_free_PRECISION( l->next_level );
      }
    } else {
      interpolation_PRECISION_dummy_free( l->next_level );
    }
    interpolation_PRECISION_free( l );
    coarse_operator_PRECISION_free( l );
  }
}


void interpolation_PRECISION_define( vector_double *V, level_struct *l, struct Thread *threading ) {

  int k, i, n = l->num_eig_vect,
      pc = 0, pi = 1, pn = n*10;
  vector_PRECISION *buffer = NULL;
  int start = threading->start_index[l->depth];
  int end   = threading->end_index[l->depth];

  if ( V == NULL ) {

    PUBLIC_MALLOC( buffer, complex_PRECISION*, 3 );
    START_MASTER(threading)
    buffer[0] = NULL;
    END_MASTER(threading)
    PUBLIC_MALLOC( buffer[0], complex_PRECISION, l->vector_size*3 );

    START_MASTER(threading)
    for( i=1; i<3; i++)
      buffer[i] = buffer[0] + l->vector_size*i;
    if ( g.print > 0 ) printf0("initial definition --- depth: %d\n", l->depth );
    if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); fflush(0); }
    if ( g.interpolation_vectors == 1 ){
      printf0("Computing eigenvectors of Q, not D! with tol %e \n", g.eigen_tol );
    }else{
      printf0("Computing eigenvectors of D, with tol %e \n", g.eigen_tol );
    }
    END_MASTER(threading)


    for ( k=0; k<n; k++ ) {
//       if ( l->depth == 0 ) {
        START_LOCKED_MASTER(threading)
        vector_PRECISION_define_random( l->is_PRECISION.test_vector[k], 0, l->inner_vector_size, l );
        END_LOCKED_MASTER(threading)
//       }

      if ( g.interpolation_vectors == 1 ){
        // multiply the test vectors with gamma5 before going into the solves
        if ( l->depth == 0 ){
          gamma5_PRECISION( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k], l, threading );
        }
        else {
          int start, end;
          compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
          SYNC_CORES(threading)
          coarse_gamma5_PRECISION( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k], start, end, l );
          SYNC_CORES(threading)
        }
      }

      // TODO : explore further error caused here by multiple OpenMP threads when GPU-enabled

      smoother_PRECISION( buffer[0], NULL, l->is_PRECISION.test_vector[k],
                          2, _NO_RES, _NO_SHIFT, l, threading );
      vector_PRECISION_copy( l->is_PRECISION.test_vector[k], buffer[0], start, end, l );

      if ( g.interpolation_vectors == 1 ){
        // multiply the test vectors with gamma5 before going into the solves
        if ( l->depth == 0 ){
          gamma5_PRECISION( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k], l, threading );
        }
        else {
          int start, end;
          compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
          SYNC_CORES(threading)
          coarse_gamma5_PRECISION( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k], start, end, l );
          SYNC_CORES(threading)
        }
      }

      smoother_PRECISION( buffer[0], NULL, l->is_PRECISION.test_vector[k],
                          g.method>=4?1:3, _NO_RES, _NO_SHIFT, l, threading );
      vector_PRECISION_copy( l->is_PRECISION.test_vector[k], buffer[0], start, end, l );

      if ( g.interpolation_vectors == 1 ){
        // multiply the test vectors with gamma5 before going into the solves
        if ( l->depth == 0 ){
          gamma5_PRECISION( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k], l, threading );
        }
        else {
          int start, end;
          compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
          SYNC_CORES(threading)
          coarse_gamma5_PRECISION( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k], start, end, l );
          SYNC_CORES(threading)
        }
      }

      smoother_PRECISION( buffer[0], NULL, l->is_PRECISION.test_vector[k],
                          g.method>=4?1:5, _NO_RES, _NO_SHIFT, l, threading );
      vector_PRECISION_copy( l->is_PRECISION.test_vector[k], buffer[0], start, end, l );

      pc += 10;
      START_MASTER(threading)
      if ( pc >= 0.2*pi*pn ) { if ( g.print > 0 ) printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); pi++; }
      END_MASTER(threading)
    }

    PUBLIC_FREE( buffer[0], complex_PRECISION, l->vector_size*3 );
    PUBLIC_FREE( buffer, complex_PRECISION*, 3 );

    for ( k=0; k<n; k++ ) {
      vector_PRECISION_real_scale( l->is_PRECISION.test_vector[k], l->is_PRECISION.test_vector[k],
                                  1.0/global_norm_PRECISION( l->is_PRECISION.test_vector[k], 0, l->inner_vector_size, l, threading ),
                                  start, end, l );
    }

    START_MASTER(threading)
    if ( g.print > 0 ) printf0("\033[0m\n");
    END_MASTER(threading)

    } else {
    for ( i=0; i<n; i++ ) {
      trans_PRECISION( l->is_PRECISION.test_vector[i], V[i], l->s_PRECISION.op.translation_table, l, threading );
    }
  }

#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
  for ( k=0; k<n; k++ ) {
    vector_PRECISION_copy( l->is_PRECISION.interpolation[k], l->is_PRECISION.test_vector[k], start, end, l );
  }
#endif



  testvector_analysis_PRECISION( l->is_PRECISION.test_vector, l, threading );

#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
  define_interpolation_PRECISION_operator( l->is_PRECISION.test_vector, l, threading );
  gram_schmidt_on_aggregates_PRECISION_vectorized( l->is_PRECISION.op, n, l, threading );
#else
  gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, n, l, threading );
  define_interpolation_PRECISION_operator( l->is_PRECISION.interpolation, l, threading );
#endif

}


void re_setup_PRECISION( level_struct *l, struct Thread *threading ) {

  if ( l->level > 0 ) {
    if ( !l->idle ) {
#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
      define_interpolation_PRECISION_operator( l->is_PRECISION.test_vector, l, threading );
      gram_schmidt_on_aggregates_PRECISION_vectorized( l->is_PRECISION.op, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_PRECISION_vectorized( l->is_PRECISION.op, l->num_eig_vect, l, threading );
      coarse_operator_PRECISION_setup_vectorized( l->is_PRECISION.op, l, threading );
      START_LOCKED_MASTER(threading)
#else
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        vector_PRECISION_copy( l->is_PRECISION.interpolation[i], l->is_PRECISION.test_vector[i],
            threading->start_index[l->depth], threading->end_index[l->depth], l );
      }
      gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, l->num_eig_vect, l, threading );
      define_interpolation_PRECISION_operator( l->is_PRECISION.interpolation, l, threading );
      START_LOCKED_MASTER(threading)
      coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
#endif
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_re_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
        }
        // else {
        coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
        //}
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_re_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
      } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
        coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
      }
      re_setup_PRECISION( l->next_level, threading );
    }
  }
//#if defined(POLYPREC) || defined(GCRODR) || defined(BLOCK_JACOBI)
#if defined(POLYPREC) || defined(GCRODR)
  else {

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_MASTER(threading)

    // this runs on level 0 only
#ifdef POLYPREC
    l->p_PRECISION.polyprec_PRECISION.update_lejas = 1;
    l->p_PRECISION.polyprec_PRECISION.preconditioner = NULL;
#endif
#ifdef GCRODR
    l->p_PRECISION.gcrodr_PRECISION.update_CU = 1;
    l->p_PRECISION.gcrodr_PRECISION.upd_ctr = 0;
    l->p_PRECISION.gcrodr_PRECISION.CU_usable = 0;
#endif
//#ifdef BLOCK_JACOBI
#if 0
    l->p_PRECISION.block_jacobi_PRECISION.local_p.polyprec_PRECISION.update_lejas = 1;
    l->p_PRECISION.block_jacobi_PRECISION.BJ_usable = 0;
#endif

    END_MASTER(threading)

    // TODO : some flags being set within this function are redundant as of
    //        the ones set above
    coarsest_level_resets_PRECISION( l, threading );

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
  }
#endif
}


void inv_iter_2lvl_extension_setup_PRECISION( int setup_iter, level_struct *l, struct Thread *threading ) {

  if ( !l->idle ) {
    vector_PRECISION buf1 = NULL;
    gmres_PRECISION_struct gmres;

    // TODO: bugfix - threading, etc

    START_LOCKED_MASTER(threading)
    MALLOC( buf1, complex_PRECISION, l->vector_size );
    fgmres_PRECISION_struct_init( &gmres );
    fgmres_PRECISION_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol,
                                   _COARSE_GMRES, _NOTHING, NULL, apply_coarse_operator_PRECISION, &gmres, l->next_level );

    if ( g.odd_even && l->next_level->level == 0 )
      gmres.v_end = l->next_level->oe_op_PRECISION.num_even_sites*l->next_level->num_lattice_site_var;
    END_LOCKED_MASTER(threading)

    for ( int k=0; k<setup_iter; k++ ) {
      int pc = 0, pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;
      START_MASTER(threading)
      printf0("depth: %d, 2lvl correction step number %d...\n", l->depth, k+1 );
      printf0("\033[0;42m\033[1;37m|"); fflush(0);
      END_MASTER(threading)
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        restrict_PRECISION( gmres.b, l->is_PRECISION.test_vector[i], l, threading );
        if ( !l->next_level->idle ) {
          if ( g.odd_even && l->next_level->level == 0 ) {
            coarse_solve_odd_even_PRECISION( &gmres, &(l->next_level->oe_op_PRECISION), l->next_level, threading );
          } else {
            fgmres_PRECISION( &gmres, l->next_level, threading );
          }
        }
        interpolate3_PRECISION( buf1, gmres.x, l, threading );
        smoother_PRECISION( buf1, NULL, l->is_PRECISION.test_vector[i], l->post_smooth_iter, _RES, _NO_SHIFT, l, threading );
        vector_PRECISION_real_scale( l->is_PRECISION.test_vector[i], buf1,
                                     1.0/global_norm_PRECISION( buf1, 0, l->inner_vector_size, l, threading ),
                                     threading->start_index[l->depth], threading->end_index[l->depth], l );
        pc += l->post_smooth_iter;
        START_MASTER(threading)
        if ( pc >= 0.2*pi*pn ) { printf0("%4d%% |", 20*pi); fflush(0); pi++; }
        END_MASTER(threading)
      }
      START_MASTER(threading)
      printf0("\033[0m\n");
      END_MASTER(threading)

#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
      define_interpolation_PRECISION_operator( l->is_PRECISION.test_vector, l, threading );
      gram_schmidt_on_aggregates_PRECISION_vectorized( l->is_PRECISION.op, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_PRECISION_vectorized( l->is_PRECISION.op, l->num_eig_vect, l, threading );
      coarse_operator_PRECISION_setup_vectorized( l->is_PRECISION.op, l, threading );
      START_LOCKED_MASTER(threading)
#else
      for ( int i=0; i<l->num_eig_vect; i++ )
        vector_PRECISION_copy( l->is_PRECISION.interpolation[i], l->is_PRECISION.test_vector[i],
            threading->start_index[l->depth], threading->end_index[l->depth], l );
      gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, l->num_eig_vect, l, threading );
      define_interpolation_PRECISION_operator( l->is_PRECISION.interpolation, l, threading );
      START_LOCKED_MASTER(threading)
      coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
#endif
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_re_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
        } else {
          coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_re_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
      } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
        coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
      }
    }

    if ( l->level > 1 )
      inv_iter_2lvl_extension_setup_PRECISION( setup_iter, l->next_level, threading );

    START_LOCKED_MASTER(threading)
    FREE( buf1, complex_PRECISION, l->vector_size );
    fgmres_PRECISION_struct_free( &gmres, l );
    END_LOCKED_MASTER(threading)
  }
}


void set_kcycle_tol_PRECISION( PRECISION tol, level_struct *l ) {

  if ( !l->idle )
    l->p_PRECISION.tol = tol;

  if ( l->level > 1 )
    set_kcycle_tol_PRECISION( tol, l->next_level );
}


void test_vector_PRECISION_update( int i, level_struct *l, struct Thread *threading ) {

  if ( l->level > 1 )
    test_vector_PRECISION_update( i, l->next_level, threading );

  if ( !l->idle && i<l->num_eig_vect ) {
#ifdef CUDA_OPT
    if( l->depth==0 ){
      vector_PRECISION_real_scale( l->is_PRECISION.test_vector[i], l->p_PRECISION.xtmp,
                                   1.0/global_norm_PRECISION( l->p_PRECISION.xtmp, 0, l->inner_vector_size, l, threading ),
                                   threading->start_index[l->depth], threading->end_index[l->depth], l );
    }
    else{
      vector_PRECISION_real_scale( l->is_PRECISION.test_vector[i], l->p_PRECISION.x,
                                   1.0/global_norm_PRECISION( l->p_PRECISION.x, 0, l->inner_vector_size, l, threading ),
                                   threading->start_index[l->depth], threading->end_index[l->depth], l );
    }
#else
    vector_PRECISION_real_scale( l->is_PRECISION.test_vector[i], l->p_PRECISION.x,
                                 1.0/global_norm_PRECISION( l->p_PRECISION.x, 0, l->inner_vector_size, l, threading ),
                                 threading->start_index[l->depth], threading->end_index[l->depth], l );
#endif
  }
}


void rayleigh_ritz_extraction_PRECISION( level_struct *l, struct Thread *threading ) {

  int m = l->num_eig_vect;

  // allocate the buffers to multiply D times the test vectors
  vector_PRECISION *buff_D_times_TVs = NULL;
  MALLOC( buff_D_times_TVs, vector_PRECISION, m );
  for ( int i=0; i < m; i++ ){
    buff_D_times_TVs[i] = NULL;
    MALLOC( buff_D_times_TVs[i], complex_PRECISION, l->inner_vector_size );
  }

  // the core of the small-subspace eigensolver
  eigslvr_PRECISION_struct eigslvr;
  eigslvr.Hc = NULL;
  eigslvr.vl = NULL;
  eigslvr.vr = NULL;
  eigslvr.w = NULL;

  // allocate memory for the compression matrix i.e. living in the small subspace
  MALLOC( eigslvr.Hc, complex_PRECISION*, m );
  eigslvr.Hc[0] = NULL;
  MALLOC( eigslvr.Hc[0], complex_PRECISION, m*m );
  for ( int i=1; i<m; i++ ) {
    eigslvr.Hc[i] = eigslvr.Hc[0] + i*m;
  }

  // allocate memory for the left and right eigenvectors
  MALLOC( eigslvr.vl, complex_PRECISION, m*m );
  MALLOC( eigslvr.vr, complex_PRECISION, m*m );
  MALLOC( eigslvr.w, complex_PRECISION, m );

  // we want to compute the right eigenvectors
  eigslvr.jobvl = 'N';
  eigslvr.jobvr = 'V';

  eigslvr.N = m;
  eigslvr.lda = m;
  eigslvr.ldvl = m;
  eigslvr.ldvr = m;
  //eigslvr.w = p->polyprec_PRECISION.h_ritz;
  //p->polyprec_PRECISION.Hc = p->polyprec_PRECISION.eigslvr.Hc;
  //p->polyprec_PRECISION.eigslvr.eigslvr_PRECISION = eigslvr_PRECISION;
  eigslvr.A = eigslvr.Hc[0];

  // apply D on the test vectors
  for ( int i=0; i < m; i++ ) {

    apply_operator_PRECISION( buff_D_times_TVs[i], l->is_PRECISION.test_vector[i], &(l->p_PRECISION), l, threading );

    // apply ganna_5 if Q is computed
    if ( g.interpolation_vectors == 1 ){
      // multiply the test vectors with gamma5 before going into the solves
      if ( l->depth == 0 ){
        gamma5_PRECISION( l->is_PRECISION.test_vector[i], l->is_PRECISION.test_vector[i], l, threading );
      }
      else {
        int start, end;
        compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
        SYNC_CORES(threading)
        coarse_gamma5_PRECISION( l->is_PRECISION.test_vector[i], l->is_PRECISION.test_vector[i], start, end, l );
        SYNC_CORES(threading)
      }
    }
  }

  // IMPORTANT : the function eigslvr_PRECISION(...) assumes LAPACK_ROW_MAJOR
  // build the small-subspace matrix
  for ( int j=0; j<m; j++ ) {
    for ( int i=0; i<m; i++ ) {
      eigslvr.Hc[i][j] = global_inner_product_PRECISION( l->is_PRECISION.test_vector[i], buff_D_times_TVs[j], 0, l->inner_vector_size, l, threading );
    }
  }

  // call the eigensolver
  eigslvr_PRECISION( &eigslvr );

  printf0("\nApproximate eigenvalues from the small subspace (info = %d):\n", eigslvr.info);
  for ( int i=0; i<m; i++ ) {
    printf0( "%f+(i)%f\n", CSPLIT(eigslvr.w[i]) );
  }
  printf0("\n");

  // finally, do the axpys to update the test vectors
  vector_PRECISION *buff_axpys = buff_D_times_TVs;
  for ( int j=0; j<m; j++ ) {
    vector_PRECISION_define( buff_axpys[j], 0, 0, l->inner_vector_size, l );
    for ( int i=0; i<m; i++ ) {
      vector_PRECISION_saxpy( buff_axpys[j], buff_axpys[j], l->is_PRECISION.test_vector[i], eigslvr.vr[j+i*m], 0, l->inner_vector_size, l );
    }
  }

  for ( int i=0; i<m; i++ ) {
    vector_PRECISION_copy( l->is_PRECISION.test_vector[i], buff_axpys[i], 0, l->inner_vector_size, l);
  }

  //printf0("RIGHT AFTER EIGENSOLVING!\n");
  //MPI_Finalize();
  //exit(0);

  // freeing memory

  for ( int i=0; i < m; i++ ){
    FREE( buff_D_times_TVs[i], complex_PRECISION, l->inner_vector_size );
  }
  FREE( buff_D_times_TVs, vector_PRECISION, m );

  FREE( eigslvr.Hc[0], complex_PRECISION, m*m );
  FREE( eigslvr.Hc, complex_PRECISION*, m );

  FREE( eigslvr.vl, complex_PRECISION, m*m );
  FREE( eigslvr.vr, complex_PRECISION, m*m );
  FREE( eigslvr.w, complex_PRECISION, m );

}



void inv_iter_inv_fcycle_PRECISION( int setup_iter, level_struct *l, struct Thread *threading ) {

  /*
    a comment on g.interpolation_vectors = 1 :

    we want to use the low modes of Q as the test vectors. This means that, instead of solving
    the block problems:

    D X = B

    we are now solving:

    Q X = B  <=>  D X = G5 B
  */
#ifdef CUDA_OPT
  if( l->depth==0 ){
    START_LOCKED_MASTER(threading)
    cuda_safe_call( cudaMallocHost( (void**)&(l->p_PRECISION.xtmp), l->inner_vector_size * sizeof(complex_PRECISION) ) );
    END_LOCKED_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
  }
#endif

  vector_PRECISION v_buf = NULL;
  complex_PRECISION *buffer = NULL;

  PUBLIC_MALLOC( buffer, complex_PRECISION, 2*l->num_eig_vect );

  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 )
    set_kcycle_tol_PRECISION( g.coarse_tol, l );
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  PUBLIC_MALLOC( v_buf, complex_PRECISION, l->vector_size );

  if ( !l->idle ) {
    //for ( int j=0; j<setup_iter; j++ ) {
    int j = 0;
    PRECISION largest_eigen_res = 1.0;
    PRECISION target_eigen_res = g.eigen_tol;

    while ( (l->depth > 0 && j < setup_iter) || (l->depth == 0 && largest_eigen_res > target_eigen_res) ) {
      int pc = 0, pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;

      START_LOCKED_MASTER(threading)
      if ( g.print > 0 ) printf0("depth: %d, bootstrap step number %d...\n", l->depth, j+1 );
      if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); if ( g.my_rank == 0 ) fflush(0); }
      END_LOCKED_MASTER(threading)

      //gram_schmidt_PRECISION( l->is_PRECISION.test_vector, buffer, 0, l->num_eig_vect, l, threading );
      //gram_schmidt_PRECISION( l->is_PRECISION.test_vector, buffer, 0, l->num_eig_vect, l, threading );

      int iters_repeat = 1;
      //if ( j > 20 ) { iters_repeat = 5; }
      for ( int ix = 0; ix < iters_repeat; ix++ ) {

        for ( int i=0; i<l->num_eig_vect; i++ ) {
          if ( g.interpolation_vectors == 1 ){
          // multiply the test vectors with gamma5 before going into the solves
          if ( l->depth == 0 ){
            gamma5_PRECISION( l->is_PRECISION.test_vector[i], l->is_PRECISION.test_vector[i], l, threading );
          }
          else {
            int start, end;
            compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
            SYNC_CORES(threading)
            coarse_gamma5_PRECISION( l->is_PRECISION.test_vector[i], l->is_PRECISION.test_vector[i], start, end, l );
            SYNC_CORES(threading)
          }
        }

#ifdef CUDA_OPT
          if( l->depth==0 ){
            vcycle_PRECISION( l->p_PRECISION.xtmp, NULL, l->is_PRECISION.test_vector[i], _NO_RES, l, threading );
          }
          else{
            vcycle_PRECISION( l->p_PRECISION.x, NULL, l->is_PRECISION.test_vector[i], _NO_RES, l, threading );
          }
#else


          if ( l->depth > 0 ) {
            vcycle_PRECISION( l->p_PRECISION.x, NULL, l->is_PRECISION.test_vector[i], _NO_RES, l, threading );
          } else {

            //int iter = 0, start = threading->start_index[l->depth], end = threading->end_index[l->depth];
            int iter = 0;
            //vector_PRECISION rhs = g.mixed_precision==2?g.p_MP.double_section.b:g.p.b;
            //vector_PRECISION sol = g.mixed_precision==2?g.p_MP.double_section.x:g.p.x;
            vector_PRECISION rhs = g.p.b;
            vector_PRECISION sol = g.p.x;
            //vector_double_copy( rhs, l->is_PRECISION.test_vector[i], start, end, l );
            trans_back_PRECISION( rhs, l->is_PRECISION.test_vector[i], l->s_PRECISION.op.translation_table, l, threading );

            int buffx = g.print;
            g.print = -1;
            iter = fgmres_PRECISION( &(g.p), l, threading );
            g.print = buffx;

            trans_PRECISION( l->p_PRECISION.x, sol, l->s_PRECISION.op.translation_table, l, threading );
            //vector_double_copy( l->p_PRECISION.x, sol, start, end, l );
            printf0("-- just did a fine-grid solve, iters = %d\n", iter);

          }

#endif
          test_vector_PRECISION_update( i, l, threading );

          pc += l->post_smooth_iter;
          START_MASTER(threading)
          if ( pc >= (int)((0.2*pi)*pn) ) { if ( g.print > 0 ) { printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); } pi++; }
          END_MASTER(threading)
        }

      }

      gram_schmidt_PRECISION( l->is_PRECISION.test_vector, buffer, 0, l->num_eig_vect, l, threading );
      //gram_schmidt_PRECISION( l->is_PRECISION.test_vector, buffer, 0, l->num_eig_vect, l, threading );

      // do the Rayleigh-Ritz extraction here!
      if ( l->depth == 0 ) {
        rayleigh_ritz_extraction_PRECISION( l, threading );
      }

      if ( l->depth == 0 ) {
        //testvector_analysis_PRECISION( l->is_PRECISION.test_vector, l, threading );
        testvector_max_residual_PRECISION( l->is_PRECISION.test_vector, &largest_eigen_res, l, threading );
      }

      START_MASTER(threading)
      if ( g.print > 0 ) printf0("\033[0m\n");
      END_MASTER(threading)

      re_setup_PRECISION( l, threading );

      if ( l->depth == 0 && l->next_level->level > 0 ) {
        inv_iter_inv_fcycle_PRECISION( MAX(1,round( ((double)(j+1)*l->next_level->setup_iter)/
        ((double)setup_iter) )), l->next_level, threading );
      }
      j++;
    }

    if ( l->depth > 0 && l->next_level->level > 0 ) {
      inv_iter_inv_fcycle_PRECISION( MAX(1, round((double)(l->next_level->setup_iter*setup_iter)/
      ((double)l->setup_iter))), l->next_level, threading );
    }
  }

  PUBLIC_FREE( v_buf, complex_PRECISION, l->vector_size );
  PUBLIC_FREE( buffer, complex_PRECISION, 2*l->num_eig_vect );

  if ( l->depth == 0 ) {
    START_LOCKED_MASTER(threading)
    set_kcycle_tol_PRECISION( g.kcycle_tol, l );
    END_LOCKED_MASTER(threading)
  }
}


void testvector_analysis_PRECISION( vector_PRECISION *test_vectors, level_struct *l, struct Thread *threading ) {
#ifdef TESTVECTOR_ANALYSIS
  START_UNTHREADED_FUNCTION(threading)
  if ( l->depth == 0 ) {

    complex_PRECISION lambda;
    PRECISION mu;
    printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );
    for ( int i=0; i<l->num_eig_vect; i++ ) {
      printf0("vector #%02d: ", i+1 );
      apply_operator_PRECISION( l->vbuf_PRECISION[3], test_vectors[i], &(l->p_PRECISION), l, no_threading );
      lambda = global_inner_product_PRECISION( test_vectors[i], l->vbuf_PRECISION[3], 0, l->inner_vector_size, l, no_threading );
      lambda /= global_inner_product_PRECISION( test_vectors[i], test_vectors[i], 0, l->inner_vector_size, l, no_threading );

      vector_PRECISION_saxpy( l->vbuf_PRECISION[1], l->vbuf_PRECISION[3], test_vectors[i], -lambda, 0, l->inner_vector_size, l );

      mu = global_norm_PRECISION( l->vbuf_PRECISION[1], 0, l->inner_vector_size, l, no_threading )/global_norm_PRECISION( test_vectors[i], 0, l->inner_vector_size, l, no_threading );
      printf0("Rayleigh quotient: %+lf%+lfi, Approximate eigenvector precision: %le\n", (double)creal(lambda), (double)cimag(lambda), (double)mu );
    }
    printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );

  }
  END_UNTHREADED_FUNCTION(threading)
#endif
}


void testvector_max_residual_PRECISION( vector_PRECISION *test_vectors, PRECISION *largest_eigen_res, level_struct *l, struct Thread *threading ) {
#ifdef TESTVECTOR_ANALYSIS
  START_UNTHREADED_FUNCTION(threading)
  PRECISION max_mu = 0.0;

  if ( l->depth == 0 ) {
    complex_PRECISION lambda;
    PRECISION mu;

    printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );
    for ( int i=0; i<l->num_eig_vect; i++ ) {
      printf0("vector #%02d: ", i+1 );

      apply_operator_PRECISION( l->vbuf_PRECISION[3], test_vectors[i],
                                &(l->p_PRECISION), l, no_threading );
      //NOTE: This function is called at finest only
      if (g.interpolation_vectors == 1) {
        gamma5_PRECISION(l->vbuf_PRECISION[3], l->vbuf_PRECISION[3], l, no_threading);
      }

      lambda = global_inner_product_PRECISION( test_vectors[i], l->vbuf_PRECISION[3],
                                               0, l->inner_vector_size, l, no_threading );
      lambda /= global_inner_product_PRECISION( test_vectors[i], test_vectors[i],
                                                0, l->inner_vector_size, l, no_threading );

      vector_PRECISION_saxpy( l->vbuf_PRECISION[1], l->vbuf_PRECISION[3], test_vectors[i],
                              -lambda, 0, l->inner_vector_size, l );

      mu = global_norm_PRECISION( l->vbuf_PRECISION[1], 0, l->inner_vector_size, l, no_threading )
         / global_norm_PRECISION( test_vectors[i], 0, l->inner_vector_size, l, no_threading );

      printf0("Rayleigh quotient: %+lf%+lfi, Approximate eigenvector precision: %le\n",
              (double)creal(lambda), (double)cimag(lambda), (double)mu );

      if (mu > max_mu) {
        max_mu = mu;
      }
    }

    printf0("Maximum residual (mu) among eigenvectors: %le\n", (double)max_mu );
    printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );
  }

  if (largest_eigen_res) {
    *largest_eigen_res = max_mu;
  }

  END_UNTHREADED_FUNCTION(threading)
#endif
}


