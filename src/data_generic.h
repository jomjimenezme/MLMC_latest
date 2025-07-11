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

#ifndef DATA_PRECISION_HEADER
  #define DATA_PRECISION_HEADER
  
  void vector_PRECISION_define( vector_PRECISION phi, complex_PRECISION value, int start, int end, level_struct *l );
  void vector_PRECISION_define_random( vector_PRECISION phi, int start, int end, level_struct *l );

  void vector_PRECISION_define_random_rademacher( vector_PRECISION phi, int start, int end, level_struct *l ); 

  void vector_PRECISION_ghg(vector_PRECISION phi, int start, int end, level_struct *l);

#endif
