#include "main.h"
#include "data_layout.h"

//TODO: We need a function that frees the memory allocated by the calls of probing, there are memory leaks

int max(int v[], int size) {
    int max = v[0];

    for (int i = 1; i < size; i++) {
        if (v[i] > max) {
            max = v[i];
        }
    }

    return max;
}

int pow_int(int base, int exp) {
    int res = 1;
    for (int i = 0; i < exp; ++i) {
        res *= base;
    }
    return res;
}

void vector_copy(int *dest, int *src, int size) {
    for (int i = 0; i < size; i++) {
        dest[i] = src[i];
    }
}

//TODO: move all the variance related functions here into a new file
void mlmc_connected_print_variances(){
  if(g.my_rank==0){
    for(int i=0; i<g.num_levels; i++){
        for(int j=0; j<g.num_levels; j++){
            int nlevs = g.num_levels;
            int idx = i*nlevs + j;
            printf("\n Variance of G_{%d,%d}(t=%d) = %f ", i,j,g.time_slice,g.variances[idx]);
        }
    }
  }
}

// Variances must be set to zero for each time-slice trace estimation
void set_probing_variances_to_zero(){
  if(g.my_rank == 0){
    if(g.trace_op_type==7){
    	for(int level = 0; level < g.num_levels*g.num_levels; level++){
            g.variances[level] = 0.0;
    	}

  }else{
	for(int level = 0; level < g.num_levels; level++){
            g.variances[level] = 0.0;
        }

    }
  }
}

//PRINT COLORS PER MPI RANK AND PER TIMESLICE
void print_colors(){

  char filename[100];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int num_processes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  sprintf(filename, "print_files/print_colors/colors_timeslice_%d_rank_%d.txt", g.time_slice, rank);
  FILE *file = fopen(filename, "w");

  fprintf(file, "\nTimeslice %d, process rank %d, t = %d, z = %d, y = %d, x = %d\n", g.time_slice, rank, g.my_coords[0], g.my_coords[1], g.my_coords[2], g.my_coords[3]);

  for(int level = 0; level < g.num_levels; level++){
    fprintf(file, "\nColors at level %d\n[", level);
    int size = g.global_lattice[level][0]*g.global_lattice[level][1]*g.global_lattice[level][2]*g.global_lattice[level][3];
    int local_size = size/num_processes;
    for(int i = 0; i < local_size; i++){
      fprintf(file, " %d", g.local_colors[level][i]);
    }
    fprintf(file, "]\n");
  }
  fclose(file); 
}

//PRINT COLORS OF ENTIRE LATTICE
void print_global_colors(){

  FILE *file = fopen("print_files/colors.txt", "w");

  for(int i = 0; i < g.num_levels; i++){
    fprintf(file, "\nColors at level %d\n [", i);

    int T = g.global_lattice[i][0];
    int Z = g.global_lattice[i][1];
    int Y = g.global_lattice[i][2];
    int X = g.global_lattice[i][3];

    int size = T*Z*Y*X;

    for(int j = 0; j < size; j++){
      fprintf(file, " %d ", g.colors[i][j]);
    }
    fprintf(file, " ]\n");
    }

    fclose(file);
}

void allocate_variances(){
//If we are doing mlmc with connected operator we have g.num_levels^2 operators
   if(g.my_rank==0){
      if(g.trace_op_type==7)
          MALLOC(g.variances, double, g.num_levels*g.num_levels);
      else
          MALLOC(g.variances, double, g.num_levels);

      set_probing_variances_to_zero();
   }
}

void setup_local_colors(){

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    // Cartisian Topology variables TODO: dims should be read from l->global_splitting
    int dims[4], periods[4], coords[4];
    MPI_Cart_get(g.comm_cart, 4, dims, periods, coords);

    g.local_colors = (int**)malloc(g.num_levels * sizeof(int*));

    for(int level = 0; level < g.num_levels; level++){

        // int size = T * Z * Y * X;
        // int local_size = size/num_processes;
        // int *current_level_colors;

        //global and local sizes in terms of cartisian grid
        int T = g.global_lattice[level][0];
        int Z = g.global_lattice[level][1];
        int Y = g.global_lattice[level][2];
        int X = g.global_lattice[level][3];

        int Nt_loc = T / dims[0];
        int Nz_loc = Z / dims[1];
        int Ny_loc = Y / dims[2];
        int Nx_loc = X / dims[3];

        const int global_size = T * Z * Y * X;
        const int local_size  = Nt_loc * Nz_loc * Ny_loc * Nx_loc;

        //if(g.my_rank == 0){
        //    MALLOC(current_level_colors, int, size);
        //    vector_copy(current_level_colors, g.colors[level], size);
        //}

        // buffer holding the FULL global colour array on every rank
        int *global_colors = NULL;
        if (g.my_rank == 0) {
            global_colors = g.colors[level];       //already allocated
        } else {
            MALLOC(global_colors, int, global_size);
        }

        // broadcast the full colour array from rank 0
        MPI_Bcast(global_colors, global_size, MPI_INT, 0, g.comm_cart);

        // allocate and fill local colour array for this level
        MALLOC(g.local_colors[level], int, local_size);

        // Initial global coordinates based on cartesian coordinates
        int t0 = coords[0] * Nt_loc;
        int z0 = coords[1] * Nz_loc;
        int y0 = coords[2] * Ny_loc;
        int x0 = coords[3] * Nx_loc;

        // loop over local linear index  idx
        int idx = 0;
        for (int lt = 0; lt < Nt_loc; lt++)
            for (int lz = 0; lz < Nz_loc; lz++)
                for (int ly = 0; ly < Ny_loc; ly++)
                    for (int lx = 0; lx < Nx_loc; lx++, idx++)
                    {
                        // global coordinates
                        int t = t0 + lt;
                        int z = z0 + lz;
                        int y = y0 + ly;
                        int x = x0 + lx;

                        // global site with t,z,y,x ordering TODO: use lex_index?
                        int gsite = ((t * Z + z) * Y + y) * X + x;
                        g.local_colors[level][idx] = global_colors[gsite];
                    }

    //ranks other than 0 allocated a temporary copy -> free it
    if (g.my_rank != 0){
        FREE(global_colors, int, global_size);
    }

        MPI_Barrier(g.comm_cart);

    }



    //if(g.my_rank==0){
    //    for(int i = 0; i < g.num_levels; i++){
	  //int T = g.global_lattice[i][0];
	  //int Z = g.global_lattice[i][1];
	  //int Y = g.global_lattice[i][2];
	  //int X = g.global_lattice[i][3];
    //      int size = T * Z * Y * X;
    //      FREE(g.colors[i], int*, size );
    //  }
    //}

    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Bcast(g.num_colors, g.num_levels, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g.num_colors, g.num_levels, MPI_INT, 0, g.comm_cart);
    MPI_Bcast(g.dilution, g.num_levels, MPI_INT, 0, g.comm_cart);
    //print_colors();
    MPI_Barrier(MPI_COMM_WORLD);

}

void get_coloring_dimension(){

  if(g.trace_op_type == 1 || g.trace_op_type == 2 || g.trace_op_type == 6 || g.trace_op_type == 7 || g.trace_op_type == 8)
    g.probing_dimension = 3;

  if(g.trace_op_type == 9 || g.trace_op_type == 10 || g.trace_op_type == 11 || g.trace_op_type == 12)
    g.probing_dimension = 4;
}
/*
void get_sigma_4D(){
  
  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;
    
    g.nc = 2;
  }
  
  if(g.coloring_distance == 2){
    g.sigma[0] = 1;
    g.sigma[1] = 2;
    g.sigma[2] = 3;
    g.sigma[3] = 4;
    
    g.nc = 10;
  }
  
  if(g.coloring_distance == 3){
    g.sigma[0] = 1;
    g.sigma[1] = 5;
    g.sigma[2] = 55;
    g.sigma[3] = 61;
    
    g.nc = 16;
  }
  
  if(g.coloring_distance == 4){
    g.sigma[0] = 1;
    g.sigma[1] = 8;
    g.sigma[2] = 12;
    g.sigma[3] = 18;
    
    g.nc = 64;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 38;
    g.sigma[1] = 1;
    g.sigma[2] = 12;
    g.sigma[3] = 16;

    g.nc = 128;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 3;
    g.sigma[1] = 20;
    g.sigma[2] = 48;
    g.sigma[3] = 50;

    g.nc = 320;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 40;
    g.sigma[1] = 32;
    g.sigma[2] = 33;
    g.sigma[3] = 61;

    g.nc = 512;
  }
  
}

void get_sigma_3D(){
  
  if(g.coloring_distance == 1){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;
    
    g.nc = 2;
  }
  
  if(g.coloring_distance == 2){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 2;
    g.sigma[3] = 3;
    
    g.nc = 8;
  }
  
  if(g.coloring_distance == 3){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 3;
    g.sigma[3] = 5;
    
    g.nc = 16;
  }
  
  if(g.coloring_distance == 4){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 6;
    g.sigma[3] = 9;
    
    g.nc = 32;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 11;
    g.sigma[3] = 27;

    g.nc = 88;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 8;
    g.sigma[3] = 44;

    g.nc = 128;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 9;
    g.sigma[3] = 33;

    g.nc = 176;
  }

  if(g.coloring_distance == 8){
    g.sigma[0] = 0;
    g.sigma[1] = 7;
    g.sigma[2] = 48;
    g.sigma[3] = 51;

    g.nc = 272;
  }

  if(g.coloring_distance == 9){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 33;
    g.sigma[3] = 45;

    g.nc = 352;
  }
  
}
*/
void dilution_check(int level){

  if(g.dilution[level] != 1 && g.dilution[level] != 2 && g.dilution[level] != 3 && g.dilution[level] != 4 && g.dilution[level] != 12){
    printf("\nError: choose a correct dilution value (1, 2, 3, 4, 12) at level %d", level);
    exit(1);
  }

  if(g.dilution[level] == 1)
    printf("\nNo dilution at level %d\n", level);

  if(g.dilution[level] == 2)
    printf("\nPartial spin dilution at level %d\n", level);

  if(g.dilution[level] == 3)
    printf("\nColor dilution at level %d\n", level);

  if(g.dilution[level] == 4)
    printf("\nComplete spin dilution at level %d\n", level);

  if(g.dilution[level] == 12)
    printf("\nSpin-Color dilution at level %d\n", level);
}
/*
void coloring_scheme(){

  MALLOC(g.num_colors, int, g.num_levels);
  MALLOC(g.dilution_ml, int, g.num_levels);

  get_coloring_dimension();

  if(g.my_rank==0){

    MALLOC(g.variances, double, g.num_levels);

    printf("\nProbing = %d - Classical probing\n", g.probing);
    printf("Coloring_distance = %d\n", g.coloring_distance);
    printf("Grids to be colored = %d\n", g.colored_grids);
    printf("Coloring dimension = %d\n", g.probing_dimension);

    if(g.probing_dimension == 3)
      get_sigma_3D();
    else
      get_sigma_4D();

    printf("sigma: %d %d %d %d\n", g.sigma[0], g.sigma[1], g.sigma[2], g.sigma[3]);
    printf("colors at the finest: %d\n", g.nc);

    dilution_check();

    double time_taken;

    double start_time = MPI_Wtime();

    g.colors = (int**)malloc(g.num_levels * sizeof(int*));

    if (g.colors == NULL)
        error0("Allocation error0\n");

    for(int level = 0; level < g.num_levels; level++){

      g.dilution_ml[level] = g.dilution;

      int T = g.global_lattice[level][0];
      int Z = g.global_lattice[level][1];
      int Y = g.global_lattice[level][2];
      int X = g.global_lattice[level][3];

      int size[4];

      size[0] = T;
      size[1] = Z;
      size[2] = Y;
      size[3] = X;

      int total_points = T * Z * Y * X;

      g.colors[level] = NULL;
      MALLOC(g.colors[level], int, total_points);
      if(level == 0){

        g.num_colors[level] = g.nc;

        // Set all colors to -1 (not assigned)
        for (int i = 0; i < total_points; i++) {
          g.colors[level][i] = -1;
        }

        // Iterate over all lattice sites
        for (int t = 0; t < T; t++) {
          for (int z = 0; z < Z; z++) {
            for (int y = 0; y < Y; y++) {
              for (int x = 0; x < X; x++) {
                int index = lex_index(t, z, y, x, size);

                // Skip if the site has already been assigned a color
                if (g.colors[level][index] != -1) {
                  continue;
                }

                //c(x) = \sum_{i = 1}^d i*x_i mod g.nc  ---> lattice dimensions labeled from 1 to d
                //int col = t + 2*z + 3*y + 4*x;
		int col = g.sigma[0]*t + g.sigma[1]*z + g.sigma[2]*y + g.sigma[3]*x;

                g.colors[level][index] = col%g.nc;

              }
            }
          }
        }

        for (int i = 0; i < total_points; i++) {
          g.colors[level][i]++;
        }


      }
      else{
        g.num_colors[level] = 1;

        for(int i = 0; i < total_points; i++)
          g.colors[level][i] = 1;
      }
    }

    double end_time = MPI_Wtime();

    time_taken = end_time - start_time;

    printf("\nTime for coloring: %f seconds\n", time_taken);
    for (int level = 0; level < g.num_levels; level++){
       printf("\n Colors at depth %d : \t %d \n", level, g.num_colors[level]);
    }


    //print_global_colors();


  }
    MPI_Barrier(MPI_COMM_WORLD);
    setup_local_colors();
}
*/

void stop_hadamard(int level){

  int size = g.global_lattice[level][0]*g.global_lattice[level][1]*g.global_lattice[level][2]*g.global_lattice[level][3];

  for(int i=0; i<size; i++){
    if(g.colors[level][i] > g.n_had[level])
      g.colors[level][i] = g.n_had[level];
  }

  g.num_colors[level] = max(g.colors[level], size);
}

void hierarchical_coloring(){

  get_coloring_dimension();

  if(g.my_rank==0){

    MALLOC(g.variances, double, g.num_levels);

    printf("\nProbing = %d - Hierarchical probing\n", g.probing);
    printf("Grids to be colored = %d\n", g.colored_grids);
    printf("Coloring dimension = %d\n", g.probing_dimension);

    double time_taken;
    double start_time = MPI_Wtime();

    g.colors = (int**)malloc(g.num_levels * sizeof(int*));

    if (g.colors == NULL)
        error0("Allocation error0\n");

    for(int level = 0; level < g.num_levels; level++){

      printf("Level %d: k = %d\n", level, g.k[level]);
      g.num_colors[level] = pow_int(2, 4*(g.k[level]-1) + 1);
      printf("Colors at level %d: %d\n", level, g.num_colors[level]);
      int Lu = pow_int(2, g.k[level]-1);
      printf("Elementary color block at level %d: %d\n", level, Lu);

      dilution_check(level);

      int T = g.global_lattice[level][0];
      int Z = g.global_lattice[level][1];
      int Y = g.global_lattice[level][2];
      int X = g.global_lattice[level][3];

      int size[4];

      size[0] = T;
      size[1] = Z;
      size[2] = Y;
      size[3] = X;

      int total_points = T * Z * Y * X;

      g.colors[level] = NULL;
      MALLOC(g.colors[level], int, total_points);
      if(level <= g.colored_grids-1){

	int *arrlc;
	MALLOC(arrlc, int, g.num_colors[level]);
        
	for(int i = 0; i < g.num_colors[level]; i++)
	  arrlc[i] = i+1; //array of all possible colors from 1 to nc


        // Set all colors to -1 (not assigned)
        for (int i = 0; i < total_points; i++) {
          g.colors[level][i] = -1;
        }

	int coords[4];
        // Iterate over all lattice sites
        for (int t = 0; t < T; t++) {
          for (int z = 0; z < Z; z++) {
            for (int y = 0; y < Y; y++) {
              for (int x = 0; x < X; x++) {
                int index = lex_index(t, z, y, x, size);

                // Skip if the site has already been assigned a color
                if (g.colors[level][index] != -1) {
                  continue;
                }

		int bx[4];
		int lx[4];
		int eo = 0;

		coords[0] = t;
		coords[1] = z;
		coords[2] = y;
		coords[3] = x;

		for(int i = 0; i < 4; i++){
		  bx[i] = coords[i]/Lu;
		  eo = eo+bx[i];
		  lx[i] = coords[i] - Lu*bx[i];
		}

		eo = eo%2;

		int idx = lx[0] + Lu*(lx[1] + Lu*(lx[2] + Lu*lx[3]));
		idx = 2*idx + eo;

		g.colors[level][index] = arrlc[idx];

              }
            }
          }
        }
        //for (int i = 0; i < total_points; i++)
          //g.colors[level][i]++;
        if(g.interrupt[level] == 1)
          stop_hadamard(level);
      }
      else{
        g.num_colors[level] = 1;

        for(int i = 0; i < total_points; i++)
          g.colors[level][i] = 1;
      }
    }

    double end_time = MPI_Wtime();

    time_taken = end_time - start_time;

    printf("\nTime for coloring: %f seconds\n", time_taken);
    //print_global_colors();

  }
    MPI_Barrier(MPI_COMM_WORLD);
    setup_local_colors();
}

void graph_coloring(){

  //if(g.probing == 1) coloring_scheme();
  if(g.probing == 2) hierarchical_coloring();

}
