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
#include "proxies/dirac_proxy_PRECISION.h"
#include "oddeven_PRECISION.h"

#if defined(POLYPREC) || defined(GMRES_POLY_EXPANSION)

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

#ifdef GMRES_POLY_EXPANSION
void apply_polyprec_jacobi_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                      operator_PRECISION_struct *op, level_struct *l,
                                      struct Thread *threading )
{
  int start, end;

  // Apply the full finest-level Dirac operator: eta = D phi
  d_plus_clover_PRECISION( eta, phi, op, l, threading );

  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

  // Apply the inverse self-coupling term: eta = C^{-1} eta = C^{-1} D phi
  diag_sc_inv_PRECISION( eta, eta, &(l->sc_op_PRECISION), l, start, end );
}
#endif

static void apply_polyprec_operator_PRECISION( vector_PRECISION output,
                                               vector_PRECISION input,
                                               gmres_PRECISION_struct *p,
                                               level_struct *l,
                                               struct Thread *threading )
{
  int start, end;

  // Apply the unrelaxed operator associated with the polynomial.
  p->polyprec_PRECISION.eval_target_operator( output, input,
                                              p->polyprec_PRECISION.target_op,
                                              l, threading );

  if ( p->shift || p->polyprec_PRECISION.omega != 1.0 ) {
    compute_core_start_end_custom(p->v_start, p->v_end, &start, &end,
                                  l, threading, l->num_lattice_site_var );

    // Include the shift before applying the relaxation factor.
    if ( p->shift )
      vector_PRECISION_saxpy( output, output, input, -p->shift,
                              start, end, l );

    // Convert the unrelaxed application into omega*A.
    if ( p->polyprec_PRECISION.omega != 1.0 )
      vector_PRECISION_scale( output, output,
                              p->polyprec_PRECISION.omega,
                              start, end, l );
  }
}

#ifdef GMRES_POLY_EXPANSION

static void polyprec_global_block_define_random_PRECISION( vector_PRECISION block, gmres_PRECISION_struct *p, level_struct *l )
{
  int rhs;
  int nrhs = p->polyprec_PRECISION.construction_nrhs;
  int vl = p->polyprec_PRECISION.syst_size;

  // Generate one random vector for every right-hand side
  for ( rhs=0; rhs<nrhs; rhs++ )
    vector_PRECISION_define_random( block + rhs*vl, p->v_start, p->v_end, l );
}


static void polyprec_global_block_scale_PRECISION( vector_PRECISION output, vector_PRECISION input, complex_PRECISION alpha, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading )
{
  int rhs, start, end;
  int nrhs = p->polyprec_PRECISION.construction_nrhs;
  int vl = p->polyprec_PRECISION.syst_size;

  compute_core_start_end_custom( p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var );

  // Scale every vector in the block
  for ( rhs=0; rhs<nrhs; rhs++ )
    vector_PRECISION_scale( output + rhs*vl, input + rhs*vl, alpha, start, end, l );
}


static void polyprec_global_block_saxpy_PRECISION( vector_PRECISION output, vector_PRECISION input, vector_PRECISION update, complex_PRECISION alpha, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading )
{
  int rhs, start, end;
  int nrhs = p->polyprec_PRECISION.construction_nrhs;
  int vl = p->polyprec_PRECISION.syst_size;

  compute_core_start_end_custom( p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var );

  // Apply the same scalar update to every vector in the block
  for ( rhs=0; rhs<nrhs; rhs++ )
    vector_PRECISION_saxpy( output + rhs*vl, input + rhs*vl, update + rhs*vl, alpha, start, end, l );
}


// Based on process_inner_product_PRECISION
static complex_PRECISION polyprec_global_block_inner_product_PRECISION( vector_PRECISION x, vector_PRECISION y, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading )
{
  int rhs;
  int nrhs = p->polyprec_PRECISION.construction_nrhs;
  int vl = p->polyprec_PRECISION.syst_size;

  complex_PRECISION local_inner_product = 0.0;

  // <X,Y>_F = sum_rhs (x_rhs)^H y_rhs
  for ( rhs=0; rhs<nrhs; rhs++ ) {

    // Get the start of the vectors
    vector_PRECISION x_rhs = x + rhs*vl;
    vector_PRECISION y_rhs = y + rhs*vl;

    // Add the threaded inner product (x_rhs)^H y_rhs on this MPI process
    local_inner_product += process_inner_product_PRECISION( x_rhs, y_rhs,  p->v_start, p->v_end, l, threading );
  }

  // Sum the local Frobenius inner products over all MPI processes
  START_MASTER(threading)
  MPI_Allreduce( &local_inner_product, &((complex_PRECISION *)threading->workspace)[0], 1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0) ? g.comm_cart : l->gs_PRECISION.level_comm );
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)

  return ((complex_PRECISION *)threading->workspace)[0];
}


static PRECISION polyprec_global_block_norm_PRECISION( vector_PRECISION x, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading )
{
  complex_PRECISION norm_squared;

  // Compute the SQUARED Frobenius norm of the block
  norm_squared = polyprec_global_block_inner_product_PRECISION( x, x, p, l, threading );

  // Return the Frobenius norm
  return (PRECISION)sqrt( creal_PRECISION(norm_squared) );
}


static void polyprec_global_block_apply_operator_PRECISION( vector_PRECISION output, vector_PRECISION input, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading )
{
  int rhs;
  int nrhs = p->polyprec_PRECISION.construction_nrhs;
  int vl = p->polyprec_PRECISION.syst_size;

  // Apply the unrelaxed target operator to every right-hand side
  for ( rhs=0; rhs<nrhs; rhs++ )
    p->polyprec_PRECISION.eval_target_operator( output + rhs*vl,input + rhs*vl, p->polyprec_PRECISION.target_op, l, threading );
}




static int polyprec_global_arnoldi_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading )
{
  int i, j;
  PRECISION norm;
  complex_PRECISION hij;

  vector_PRECISION rhs = p->polyprec_PRECISION.global_rhs;
  vector_PRECISION *V = p->polyprec_PRECISION.global_V;
  vector_PRECISION w = p->polyprec_PRECISION.global_w;
  complex_PRECISION **H = p->polyprec_PRECISION.Hc;

  // Compute ||B||_F
  norm = polyprec_global_block_norm_PRECISION( rhs, p, l, threading );

  // V_0 = B / ||B||_F
  polyprec_global_block_scale_PRECISION( V[0], rhs, 1.0/norm, p, l, threading );

  for ( j=0; j<p->polyprec_PRECISION.d_poly; j++ ) {

    // W = Ahat V_j
    polyprec_global_block_apply_operator_PRECISION( w, V[j], p, l, threading );

    // Orthogonalize W against V_0,...,V_j
    for ( i=0; i<=j; i++ ) {

      // h_{i,j} = <V_i,W>_F
      hij = polyprec_global_block_inner_product_PRECISION( V[i], w, p, l, threading );

      // Store the Arnoldi coefficient
      START_MASTER(threading)
      H[j][i] = hij;
      END_MASTER(threading)

      SYNC_MASTER_TO_ALL(threading)

      // W = W - h_{i,j} V_i
      polyprec_global_block_saxpy_PRECISION( w, w, V[i], -H[j][i], p, l, threading );
    }

    // h_{j+1,j} = ||W||_F
    norm = polyprec_global_block_norm_PRECISION( w, p, l, threading );

    // Store the norm of the new Arnoldi vector
    START_MASTER(threading)
    H[j][j+1] = norm;
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading)

    // V_{j+1} = W / h_{j+1,j}
    polyprec_global_block_scale_PRECISION( V[j+1], w, 1.0/H[j][j+1], p, l, threading );
  }

  return 0;
}

#endif

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

static void finalize_polyprec_roots_PRECISION( gmres_PRECISION_struct *p )
{
  int i;

  // Compute the harmonic Ritz values from the Arnoldi Hessenberg matrix
  harmonic_ritz_PRECISION( p );

  // Scale the harmonic Ritz values for the relaxed operator omega*A
  if ( p->polyprec_PRECISION.omega != 1.0 ) {
    for ( i=0; i<p->polyprec_PRECISION.d_poly; i++ )
      p->polyprec_PRECISION.h_ritz[i] *=
        p->polyprec_PRECISION.omega;
  }

  // Order the roots
  leja_ordering_PRECISION( p );
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
  finalize_polyprec_roots_PRECISION( p );
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



#ifdef POLYPREC

int re_construct_lejas_PRECISION( level_struct *l, struct Thread *threading ) {

  //printf0("UPDATED LEJAS\n");

  return update_lejas_PRECISION(&(l->p_PRECISION), l, threading);

}

#endif

#ifdef GMRES_POLY_EXPANSION
int construct_fine_polyprec_PRECISION( gmres_PRECISION_struct *p,
                                       level_struct *l,
                                       struct Thread *threading )
{
  // Only the Jacobi splitting is implemented so far
  if ( p->polyprec_PRECISION.splitting != _POLYPREC_JACOBI )
    error0("POLYPREC: the selected finest-level splitting is not implemented yet.\n");

  START_LOCKED_MASTER(threading)

  // Construct the factors required to apply C^{-1} (stored in op->clover)
  selfcoupling_setup_PRECISION( &(g.op_double), l );

  // D (g.op_PRECISION) gives the operator data for the Jacobi application C^{-1}D
  p->polyprec_PRECISION.target_op = &(g.op_PRECISION);

  // Select the function that applies the Jacobi preconditioned operator
  p->polyprec_PRECISION.eval_target_operator = apply_polyprec_jacobi_PRECISION;

  // The roots must be constructed for the current operator and splitting!!
  p->polyprec_PRECISION.update_lejas = 1;

  END_LOCKED_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // Construct the Leja roots for the assigned finest target operator
  int polyprec_status = update_lejas_PRECISION( p, l, threading );
  // Abort because the finest polynomial is required by the expansion!!
  if ( polyprec_status != 1 )
    error0("POLYPREC: finest-level polynomial construction failed.\n");

  // Report SUCCESFUL polynomial construction
  return polyprec_status;
}
#endif

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

// Apply the complete inverse polynomial omega q_{d-1}(A_omega).
void apply_polyprec_inverse_core_PRECISION( vector_PRECISION phi, vector_PRECISION eta,
                                            gmres_PRECISION_struct *p, level_struct *l,
                                            struct Thread *threading )
{
  int start, end;
  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

  // Apply q_{d-1}(A_omega) to eta
  apply_polyprec_core_PRECISION( phi, eta, p, l, threading );

  // Apply the relaxation factor
  if ( p->polyprec_PRECISION.omega != 1.0 )
    vector_PRECISION_scale( phi, phi, p->polyprec_PRECISION.omega, start, end, l );
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


#ifdef POLYPREC

void apply_polyprec_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                               int res, level_struct *l, struct Thread *threading )
{
  apply_polyprec_core_PRECISION( phi, eta, &(l->p_PRECISION), l, threading );
}

#endif

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
