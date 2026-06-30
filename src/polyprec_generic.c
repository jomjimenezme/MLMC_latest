/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
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

#ifdef POLYPREC


/*-----------------------------------------------*/

void print_matrix_PRECISION(complex_PRECISION* A, int mv, int mh )
{
  int i,j;

  // printf("\n\n");
  // for (i=0; i < mv; i++)
  // {
  //     for(j=0; j < mh; j++)
  //     {
  //             fprintf(stdout, "%6.6f +i%6.6f\t", creal(A[i*mh + j]), cimag(A[i*mh+j]));
  //     }
  //     fprintf(stdout, "\n");
  // }
  // printf("--\n");

  printf("\n\n");
  for (i=0; i < mh; i++)
  {
    for(j=0; j < mv; j++)
    {
      fprintf(stdout, "%6.6f +i%6.6f\t", creal(A[j*mh + i]), cimag(A[j*mh+i]));
    }
    fprintf(stdout, "\n");
  }
  printf("--\n");
  printf("\n\n");
}


void print_vector_PRECISION( char* desc, vector_PRECISION w, int n)
{
  int j;
  printf0( "\n %s\n", desc );
  for( j = 0; j < n; j++ ) printf0( " (%6.6f,%6.6f)", creal(w[j]), cimag(w[j]) );
  printf0( "\n" );
}

static void apply_polyprec_operator_PRECISION( vector_PRECISION output,
                                               vector_PRECISION input,
                                               gmres_PRECISION_struct *p,
                                               level_struct *l,
                                               struct Thread *threading )
{
  // Apply the operator for which the polynomial was constructed.
  p->polyprec_PRECISION.eval_target_operator( output, input,
                                              p->polyprec_PRECISION.target_op,
                                              l, threading );

  if ( p->shift ) {
    int start, end;
    compute_core_start_end_custom(p->v_start, p->v_end, &start, &end,
                                  l, threading, l->num_lattice_site_var );
    vector_PRECISION_saxpy( output, output, input, -p->shift, start, end, l );
  }
}

void harmonic_ritz_PRECISION( gmres_PRECISION_struct *p )
{
  int i, j, d;
  complex_PRECISION h_dd;

  d = p->polyprec_PRECISION.d_poly;
  h_dd = p->polyprec_PRECISION.Hc[d-1][d];
  memset(p->polyprec_PRECISION.dirctslvr.b, 0.0, sizeof(complex_PRECISION)*(d-1));
  p->polyprec_PRECISION.dirctslvr.b[d-1] = 1.;

  for (i=0; i<d; i++)
    for (j=0; j<d; j++)
      p->polyprec_PRECISION.Hcc[i*d + j ] = conj(p->polyprec_PRECISION.Hc[j][i]);

  p->polyprec_PRECISION.dirctslvr.dirctslvr_PRECISION(&p->polyprec_PRECISION.dirctslvr);

  for (i=0; i<d; i++)
    p->polyprec_PRECISION.Hc[d-1][i] += h_dd*h_dd*p->polyprec_PRECISION.dirctslvr.x[i];
    
  p->polyprec_PRECISION.eigslvr.eigslvr_PRECISION(&p->polyprec_PRECISION.eigslvr);
}



/*-----------------------------------------------*/



void leja_ordering_PRECISION( gmres_PRECISION_struct *p )
{

  int i, j, ii, d_poly;
  int max_j, exchange_cols;
  complex_PRECISION tmp, leja;

  complex_PRECISION** L;
  complex_PRECISION* col_prods;

  d_poly = p->polyprec_PRECISION.d_poly;
  L = p->polyprec_PRECISION.L;
  col_prods = p->polyprec_PRECISION.col_prods;

  // Create a matrix made of n+1 rows, each row is x (all rows equal).
  for (i=0; i<d_poly+1; i++ )
    memcpy( L[i], p->polyprec_PRECISION.h_ritz, sizeof(complex_PRECISION)*(d_poly) );

  leja = 0; 

  for (i=0; i < d_poly-1; i++)
  {
    for (j=i; j<d_poly; j++ ) 
      L[i][j] = cabs( L[i][j] - leja );

    for (j = i; j < d_poly; j++)
    {
      col_prods[j] = 1.;
      for (ii = 0; ii <= i; ii++)
        col_prods[j] *= L[ii][j];
    }
        
    exchange_cols = 0;
    max_j = i;
    for (j=i+1; j<d_poly; j++ )
    {
      if ( creal(col_prods[j]) > creal(col_prods[max_j]) )
      {
        max_j = j; 
        exchange_cols = 1;
      }
    }
        
    if (exchange_cols)
    {
      for (ii=0; ii<d_poly+1; ii++ )
      {
        tmp = L[ii][i];
        L[ii][i] = L[ii][max_j];
        L[ii][max_j] = tmp;
      } 
    }

    leja = L[d_poly][i];

  }

  memcpy( p->polyprec_PRECISION.lejas, p->polyprec_PRECISION.L[d_poly], sizeof(complex_PRECISION)*(d_poly) );
}



int update_lejas_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading )
{
  int start, end;
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  vector_PRECISION random_rhs, buff0;
  random_rhs = p->polyprec_PRECISION.random_rhs;
  PRECISION buff3, buff5;
  vector_PRECISION buff4;

  int buff1, buff2, buff_initial_guess_zero;
  int fgmres_itersx;
  void (*buff_preconditioner)();

  buff0 = p->b;
  buff2 = p->num_restart;
  buff1 = p->restart_length;
  buff3 = p->tol;
  buff4 = p->x;
  buff5 = g.coarse_tol;

  // For polynomial expansion only
  operator_PRECISION_struct *buff_op;
  void (*buff_eval_operator)(vector_PRECISION, vector_PRECISION,
                             operator_PRECISION_struct *,
                             struct level_struct *, struct Thread *);

  buff_op = p->op;
  buff_eval_operator = p->eval_operator;

  // Save the initial guess and preconditioner of the GMRES workspace (p)
  // Polynomial construction temporarily changes these GMRES settings,
  // so they must be restored before this function returns.
  buff_initial_guess_zero = p->initial_guess_zero;
  buff_preconditioner = p->preconditioner;

  if ( p->polyprec_PRECISION.d_poly > buff1 ) {
  START_MASTER(threading)
  error0(
      "POLYPREC: polynomial degree %d exceeds the GMRES "
      "restart length %d used to allocate the Arnoldi workspace.\n",
      p->polyprec_PRECISION.d_poly,
      buff1
  );
  END_MASTER(threading)

  return -2;
}

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  START_MASTER(threading)
  p->b = random_rhs;
  p->num_restart = 1;
  p->restart_length = p->polyprec_PRECISION.d_poly;
  p->preconditioner = NULL;
  // Do not use the temporary solution vector as an initial guess.
  p->initial_guess_zero = 1;
  p->tol = 1E-20;
  if ( l->level == 0 )
    g.coarse_tol = 1E-20;
  p->x = p->polyprec_PRECISION.xtmp;
  // l->dup_H = 1;  (checks if Arnoldi must copy H)
  p->polyprec_PRECISION.capture_H = 1;

  // Use the operator associated with the polynomial
  p->op = p->polyprec_PRECISION.target_op;
  p->eval_operator = p->polyprec_PRECISION.eval_target_operator;

  vector_PRECISION_define_random( random_rhs, p->v_start, p->v_end, l );
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  fgmres_itersx = fgmres_PRECISION(p, l, threading);

  //printf0( "FROM WITHIN POLYPREC SETUP : %d, d POLY = %d\n",fgmres_itersx,p->polyprec_PRECISION.d_poly );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  START_MASTER(threading)
  p->polyprec_PRECISION.capture_H = 0;
  p->b = buff0;
  p->num_restart = buff2;
  p->restart_length = buff1;
  p->tol = buff3;
  if ( l->level == 0 )
    g.coarse_tol = buff5;
  p->x = buff4;
  // Restore the original GMRES state.
  p->initial_guess_zero = buff_initial_guess_zero;

  p->preconditioner = buff_preconditioner;
  p->op = buff_op;

  // Restore the operator used by the original GMRES workspace
  p->eval_operator = buff_eval_operator;
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading);

  if ( fgmres_itersx == p->polyprec_PRECISION.d_poly ) {
    START_MASTER(threading)
    p->polyprec_PRECISION.preconditioner = p->polyprec_PRECISION.preconditioner_bare;
    /*
     * Before: l->p_PRECISION.polyprec_PRECISION.update_lejas = 0;
     * but we will need to access p = &g.p; for the polynomial at
     * the finest
     */
    p->polyprec_PRECISION.update_lejas = 0;
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading);

  } else { return -1; }

  START_MASTER(threading)
  harmonic_ritz_PRECISION(p);
  leja_ordering_PRECISION(p);
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

#ifdef POLYPREC_CHECK
  PRECISION polyprec_error = check_polyprec_identity_PRECISION(p, l, threading);

  START_MASTER(threading)
  printf0("POLYPREC: polynomial identity error, p_d(A) eta = eta - A q_{d-1}(A) eta. Error: %le\n",
          polyprec_error);
  END_MASTER(threading)
#endif

  return 1;
}



int re_construct_lejas_PRECISION( level_struct *l, struct Thread *threading ) {

  //printf0("UPDATED LEJAS\n");

  return update_lejas_PRECISION(&(l->p_PRECISION), l, threading);

}


void apply_polyprec_core_PRECISION( vector_PRECISION phi, vector_PRECISION eta,
                                    gmres_PRECISION_struct *p, level_struct *l,
                                    struct Thread *threading )
{

  //printf0("APPLYING POLYNOMIAL\n");

  int i, start, end;

  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  // Polynomial degree
  int d_poly = p->polyprec_PRECISION.d_poly;
  // Accumulates the value of q_{d-1}(A) eta.
  vector_PRECISION accum_prod = p->polyprec_PRECISION.accum_prod;
  // stores the current product of residual factors ( I - A/theta_j )
  vector_PRECISION product = p->polyprec_PRECISION.product;
  // stores A times the current product
  vector_PRECISION temp = p->polyprec_PRECISION.temp;
  vector_PRECISION lejas = p->polyprec_PRECISION.lejas;

  // Initialize the first product with eta
  vector_PRECISION_copy( product, eta, start, end, l );
  // Initialize the accumulated polynomial with zero:
  vector_PRECISION_define(accum_prod, 0.0, start, end, l);
  // accum_prod = eta/theta_0.
  vector_PRECISION_saxpy(accum_prod, accum_prod, product, 1./lejas[0], start, end, l);

  //Accumulate the remaining terms i = 1,...,d-1.
  for (i = 1; i < d_poly; i++)
  {
#ifdef PERS_COMMS
    g.pers_comms_id2 = p->restart_length + g.pers_comms_nrZxs;
    g.use_pers_comms1 = 1;
#endif

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    // temp = A product
    //apply_operator_PRECISION(temp, product, p, l, threading);
    apply_polyprec_operator_PRECISION(temp, product, p, l, threading);
#ifdef PERS_COMMS
    g.pers_comms_id2 = -1;
    g.use_pers_comms1 = 0;
#endif

    // Get next residual factor: (I - A/theta_{i-1}) product
    vector_PRECISION_saxpy(product, product, temp, -1./lejas[i-1], start, end, l);
    /// Add the next term of the inverse-approximation polynomial:  accum_prod + product/theta_i
    vector_PRECISION_saxpy(accum_prod, accum_prod, product, 1./lejas[i], start, end, l);
  }

  vector_PRECISION_copy( phi, accum_prod, start, end, l );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
}

void apply_polyprec_residual_core_PRECISION( vector_PRECISION phi, vector_PRECISION eta,
                                             gmres_PRECISION_struct *p, level_struct *l,
                                             struct Thread *threading )
{
  // Evaluate the residual polynomial p_d(A) eta.
  int i, start, end;

  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  int d_poly = p->polyprec_PRECISION.d_poly;
  vector_PRECISION product = p->polyprec_PRECISION.product;
  vector_PRECISION temp = p->polyprec_PRECISION.temp;
  vector_PRECISION lejas = p->polyprec_PRECISION.lejas;

  vector_PRECISION_copy( product, eta, start, end, l );

  for (i = 0; i < d_poly; i++)
  {
#ifdef PERS_COMMS
    g.pers_comms_id2 = p->restart_length + g.pers_comms_nrZxs;
    g.use_pers_comms1 = 1;
#endif
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    //apply_operator_PRECISION(temp, product, p, l, threading);
    apply_polyprec_operator_PRECISION(temp, product, p, l, threading);
#ifdef PERS_COMMS
    g.pers_comms_id2 = -1;
    g.use_pers_comms1 = 0;
#endif

    vector_PRECISION_saxpy(product, product, temp, -1./lejas[i], start, end, l);
  }

  vector_PRECISION_copy( phi, product, start, end, l );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
}


void apply_polyprec_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                               int res, level_struct *l, struct Thread *threading )
{
  apply_polyprec_core_PRECISION( phi, eta, &(l->p_PRECISION), l, threading );
}

#ifdef POLYPREC_CHECK
PRECISION check_polyprec_identity_PRECISION( gmres_PRECISION_struct *p, level_struct *l,
                                             struct Thread *threading )
  // Check if p_d(A) eta = eta - A q_{d-1}(A) eta
{
  int start, end;
  PRECISION norm_eta, norm_diff;

  vector_PRECISION eta = p->polyprec_PRECISION.random_rhs;
  vector_PRECISION p_eta = p->polyprec_PRECISION.xtmp;
  vector_PRECISION q_eta = p->polyprec_PRECISION.accum_prod;
  vector_PRECISION check = p->polyprec_PRECISION.product;
  vector_PRECISION temp = p->polyprec_PRECISION.temp;

  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  apply_polyprec_residual_core_PRECISION( p_eta, eta, p, l, threading );
  apply_polyprec_core_PRECISION( q_eta, eta, p, l, threading );

  //apply_operator_PRECISION(temp, q_eta, p, l, threading);
  apply_polyprec_operator_PRECISION(temp, q_eta, p, l, threading);

  vector_PRECISION_copy( check, eta, start, end, l );
  vector_PRECISION_saxpy(check, check, temp, -1.0, start, end, l);
  vector_PRECISION_saxpy(check, check, p_eta, -1.0, start, end, l);

  norm_eta = global_norm_PRECISION( eta, p->v_start, p->v_end, l, threading );
  norm_diff = global_norm_PRECISION( check, p->v_start, p->v_end, l, threading );

  return norm_diff/norm_eta;
}
#endif
#endif
