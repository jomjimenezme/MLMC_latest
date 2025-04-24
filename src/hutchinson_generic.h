#ifndef HUTCHINSON_PRECISION_HEADER
  #define HUTCHINSON_PRECISION_HEADER

  struct Thread;

  void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );
  void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );
  int apply_solver_PRECISION( level_struct* l, struct Thread *threading );
   
  gmres_PRECISION_struct* get_p_struct_PRECISION( level_struct* l );

  complex_PRECISION hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );

  struct sample hutchinson_blind_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading );

  void hutchinson_diver_PRECISION_init( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading );
  complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );


#endif
