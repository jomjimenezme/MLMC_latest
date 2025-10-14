#ifndef GLOBAL_STRUCT_H
#define GLOBAL_STRUCT_H

#ifndef IMPORT_FROM_EXTERN_C
#include <mpi.h>
#endif

#include <stdio.h>
#include "global_defs.h"
#include "algorithm_structs_double.h"
#include "algorithm_structs_float.h"
#include "algorithm_structs.h"
#include "var_table.h"
#include "dd_alpha_amg_parameters.h"
#include "dd_alpha_amg_setup_status.h"

typedef struct global_struct
{

    FILE *logfile;

    gmres_double_struct p;
    gmres_MP_struct p_MP;
    operator_double_struct op_double;
    operator_float_struct op_float;

    // communication
    MPI_Comm comm_cart;
    MPI_Group global_comm_group;
    MPI_Request sreqs[8], rreqs[8];
    int num_processes, my_rank, my_coords[4], two_cnfgs, tv_io_single_file, num_openmp_processes;
    // string buffers
    char in[STRINGLENGTH], in_clov[STRINGLENGTH], source_list[STRINGLENGTH], tv_io_file_name[STRINGLENGTH];
    // geometry, method parameters
    int num_levels, num_desired_levels, process_grid[4], in_format,
        **global_lattice, **local_lattice, **block_lattice,
        *post_smooth_iter, *block_iter, *setup_iter, *ncycle,
        method, odd_even, anti_pbc, rhs, propagator_coords[4],
        interpolation, randomize, *num_eig_vect, num_coarse_eig_vect, kcycle, mixed_precision,
        restart, max_restart, kcycle_restart, kcycle_max_restart, coarse_iter, coarse_restart,
        *trace_max_iters, *trace_min_iters, time_slice, trace_op_type, time_slice_inner_connected;
    double tol, coarse_tol, kcycle_tol, csw, rho, *relax_fac;

    // Improved setup 1: yes, 0:no
    int default_setup;
    // Eigentol in improved Setup
    double eigen_tol;
    // test vectors: 0: from D,  1; from Q
    int interpolation_vectors;
    // store the test vectors 1: yes, 0:no
    int write_tv;

    // profiling, analysis, output
    int coarse_iter_count, iter_count, iterator, print, conf_flag, setup_flag, in_setup;
    double coarse_time, prec_time, *output_table[8], cur_storage, max_storage, total_time,
        plaq_hopp, plaq_clov, norm_res, plaq, setup_m0, solve_m0, bicgstab_tol;

#ifdef CUDA_OPT
    double cur_gpu_storage, max_gpu_storage;

    // some CUDA-specific values
    int warp_size;
    int num_devices;
    int device_id;

    int *CUDA_threads_per_CUDA_block_type1;
    int *CUDA_threads_per_lattice_site_type1;
    int *CUDA_threads_per_CUDA_block_type2;
    int *CUDA_threads_per_lattice_site_type2;
#endif

    // index functions for external usage
    int (*conf_index_fct)(), (*vector_index_fct)();
    int *odd_even_table;

    // bc: 0 dirichlet, 1 periodic, 2 anti-periodic
    int bc;

    complex_double **gamma, g5D_shift;
    var_table vt;

    int on_solve;
    int nr_threads;

    struct dd_alpha_amg_parameters amg_params;
    struct dd_alpha_amg_setup_status mg_setup_status;
    double mass_for_next_solve;

#ifdef CUDA_OPT
    int oddeven_copy_nt_2_gpu;
#endif

#ifdef TM_COARSEST
    double mu_coarsest;
#endif

#ifdef GCRODR
    int gcrodr_k, gcrodr_k_setup, gcrodr_k_solve;
    int gcrodr_upd_itrs_solve;
    int gcrodr_upd_itrs_setup;
    int gcrodr_calling_from_setup;
#endif

#ifdef POLYPREC
    int polyprec_d, polyprec_d_setup, polyprec_d_solve;
#endif

//#ifdef BLOCK_JACOBI
#if 0
    int local_polyprec_d;
#endif

//#ifdef BLOCK_JACOBI
#if 0
    double bj_time;
#endif
#ifdef GCRODR
    double gcrodr_LSP_time, gcrodr_buildAB_time, gcrodr_buildCU_time;
#endif

#ifdef RICHARDSON_SMOOTHER
    int smoother_richardson_BPI_iters, richardson_sub_degree;
    double richardson_factor;
#endif

    double coarsest_time;
    double matmul_time;

    double avg_b1;
    double avg_b2;
    double avg_crst;


    int probing; //contains information on whether probing is performed or not
    int probing_dimension; // 4D or 3D coloring
    int coloring_distance;
    int coloring_method;

    int **colors; //colors of the lattice
    int **local_colors; //colors of every MPI process
    int *num_colors; //number of colors at every level
    int coloring_count;
    
    int dilution;
    int *dilution_ml;
    int dilution_count;

    int nc;
    int sigma[4];

    double *variances; //variance of the estimator at every level

} global_struct;

extern global_struct g;

#endif // GLOBAL_STRUCT_H
