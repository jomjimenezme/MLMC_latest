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

#ifndef POLYPREC_PRECISION_HEADER
  #define POLYPREC_PRECISION_HEADER

#ifdef POLYPREC

  void apply_polyprec_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                                 int res, level_struct *l, struct Thread *threading );

  int re_construct_lejas_PRECISION( level_struct *l, struct Thread *threading );

#endif

#ifdef GMRES_POLY_EXPANSION

  void apply_polyprec_jacobi_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                        operator_PRECISION_struct *op, level_struct *l,
                                        struct Thread *threading );



  int construct_fine_polyprec_PRECISION( gmres_PRECISION_struct *p,
                                         level_struct *l,
                                         struct Thread *threading );

#endif

  void apply_polyprec_residual_core_PRECISION( vector_PRECISION phi, vector_PRECISION eta,
                                               gmres_PRECISION_struct *p, level_struct *l,
                                               struct Thread *threading );

  void apply_polyprec_core_PRECISION( vector_PRECISION phi, vector_PRECISION eta,
                                      gmres_PRECISION_struct *p, level_struct *l,
                                      struct Thread *threading );

  // Apply the complete inverse polynomial omega q_{d-1}(A_omega).
  void apply_polyprec_inverse_core_PRECISION( vector_PRECISION phi, vector_PRECISION eta,
                                              gmres_PRECISION_struct *p, level_struct *l,
                                              struct Thread *threading );

  #ifdef POLYPREC_CHECK
  PRECISION check_polyprec_identity_PRECISION( gmres_PRECISION_struct *p, level_struct *l,
                                               struct Thread *threading );
  #endif

#endif
