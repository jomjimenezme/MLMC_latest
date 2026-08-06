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


        //g.local_colors[level] = NULL;
        //MALLOC(g.local_colors[level], int, local_size);

        //int *current_level_local_colors;
        //MALLOC(current_level_local_colors, int, local_size);

        //MPI_Barrier(MPI_COMM_WORLD);

        //MPI_Scatter(current_level_colors, local_size, MPI_INT,
        //        current_level_local_colors, local_size, MPI_INT,
        //        0, MPI_COMM_WORLD);

        //vector_copy(g.local_colors[level], current_level_local_colors, local_size);

        //FREE(current_level_colors, int, size);
        //FREE(current_level_local_colors, int, local_size);

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


    // TODO: make it a function or remove this block
 /*   // Debug print: global colour array on rank 0
        if (g.my_rank == 0) {
            printf("[Rank %d] Global color array at level %d:\n", g.my_rank, level);
            for (int i = 0; i < global_size; i++)
                printf("%d ", global_colors[i]);
            printf("\n");
        }
*/
        MPI_Barrier(g.comm_cart);

  /*  // Debug print: local colour arrays
        for (int r = 0; r < num_processes; r++) {
            if (g.my_rank == r) {
                printf("[Rank %d] Local color array at level %d:\n", g.my_rank, level);
                for (int i = 0; i < local_size; i++)
                    printf("%d ", g.local_colors[level][i]);
                printf("\n");
                fflush(stdout);
            }
            MPI_Barrier(g.comm_cart);
        }
*/


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
    MPI_Bcast(g.dilution_ml, g.num_levels, MPI_INT, 0, g.comm_cart);
    //print_colors();
    MPI_Barrier(MPI_COMM_WORLD);

}

void generate_neighbors(int t, int z, int y, int x, int **neighbors, int *num_neighbors, int size[4]) {
    int T = size[0];
    int Z = size[1];
    int Y = size[2];
    int X = size[3];

    int max_neighbors = 1; // Includes central point
    for (int delta = 1; delta <= g.coloring_distance; delta++) {
        max_neighbors += 8 * delta * delta * delta; // Estimate maximum number
    }
    *neighbors = (int *)malloc(max_neighbors * sizeof(int));
    *num_neighbors = 0;

    // Iterate over all displacements combination
    for (int dt = -g.coloring_distance; dt <= g.coloring_distance; dt++) {
        for (int dz = -g.coloring_distance; dz <= g.coloring_distance; dz++) {
            for (int dy = -g.coloring_distance; dy <= g.coloring_distance; dy++) {
                for (int dx = -g.coloring_distance; dx <= g.coloring_distance; dx++) {

                    //Verify that the distance is the right one
                    if (abs(dt) + abs(dz) + abs(dy) + abs(dx) > g.coloring_distance) {
                        continue;
                    }

                    // Skip central point (0, 0, 0, 0)
                    if (dt == 0 && dz == 0 && dy == 0 && dx == 0) {
                        continue;
                    }

                    // Compute coordinates of the neighbor
                    int t_neighbor = (t + dt + T) % T;
                    int z_neighbor = (z + dz + Z) % Z;
                    int y_neighbor = (y + dy + Y) % Y;
                    int x_neighbor = (x + dx + X) % X;

                    // Compute the lexicographic index of the neighbor
                    (*neighbors)[*num_neighbors] = lex_index(t_neighbor, z_neighbor, y_neighbor, x_neighbor, size);
                    (*num_neighbors)++;
                }
            }
        }
    }
}

/*
// GREEDY COLORING
void graph_coloring() {

    MALLOC(g.num_colors, int, g.num_levels);

    if(g.my_rank == 0){

    printf("\nProbing = %d\n", g.probing);
    printf("Coloring_distance = %d\n", g.coloring_distance);
    printf("Coloring_method = %d\n", g.coloring_method);

    double time_taken;

    double start_time = MPI_Wtime();

    g.colors = (int**)malloc(g.num_levels * sizeof(int*));

    if (g.colors == NULL)
        error0("Allocation error0\n");

    for(int level = 0; level < g.num_levels; level++){

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

    if(level<g.coloring_method){

    // Set all colors to 0 (not assigned)
    for (int i = 0; i < total_points; i++) {
        g.colors[level][i] = 0;
    }

    // Iterate over all lattice sites
    for (int t = 0; t < T; t++) {
        if (g.probing_dimension == 3 && t != g.time_slice) continue;  //TODO: This should be done for 3D traces ONLY!!
        for (int z = 0; z < Z; z++) {
            for (int y = 0; y < Y; y++) {
                for (int x = 0; x < X; x++) {
                    int index = lex_index(t, z, y, x, size);

                    // Skip if the site has already been assigned a color
                    if (g.colors[level][index] != 0) {
                        continue;
                    }

                    // Find neighbors within distance d
                    int *neighbors;
                    int num_neighbors;
                    generate_neighbors(t, z, y, x, &neighbors, &num_neighbors, size);

                    int used_colors[256] = {0};
                    for (int i = 0; i < num_neighbors; i++) {
                        int neighbor_index = neighbors[i];
                        if (g.colors[level][neighbor_index] != 0) {
                            used_colors[g.colors[level][neighbor_index] - 1] = 1;
                        }
                    }

                    // Assign the minimum possible color
                    int color = 1;
                    while (used_colors[color - 1]) {
                        color++;
                    }
                    g.colors[level][index] = color;

                    // Free memory of neighbors
                    free(neighbors);
                }
            }
        }
    }
    }
    else{
        for(int i = 0; i < total_points; i++)
            g.colors[level][i] = 1;
        }

    g.num_colors[level] = max(g.colors[level], total_points);

    }

    double end_time = MPI_Wtime();

    time_taken = end_time - start_time;

    printf("\nTime for coloring: %f seconds\n", time_taken);
    for (int level = 0; level < g.num_levels; level++){
       printf("\n Colors at depth %d : \t %d \n", level, g.num_colors[level]);
    }

    /*
    FILE *file = fopen("print_files/colors.txt", "w");

    for(int i = 0; i < g.num_levels; i++){
        fprintf(file, "\nColors at level %d\n [", i+1);

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
    */ /*
    }

    MPI_Barrier(MPI_COMM_WORLD);
    setup_local_colors();
}

*/

void get_sigma_4D(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.nc = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 4;
    g.sigma[1] = 1;
    g.sigma[2] = 9;
    g.sigma[3] = 10;

    g.nc = 15;
  }

  if(g.coloring_distance == 3){
    g.sigma[0] = 19;
    g.sigma[1] = 7;
    g.sigma[2] = 21;
    g.sigma[3] = 15;

    g.nc = 24;
  }

  if(g.coloring_distance == 4){
    g.sigma[0] = 33;
    g.sigma[1] = 45;
    g.sigma[2] = 16;
    g.sigma[3] = 7;

    g.nc = 88;
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
    g.sigma[2] = 4;
    g.sigma[3] = 10;

    g.nc = 32;
  }

  if(g.coloring_distance == 5){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 5;
    g.sigma[3] = 13;

    g.nc = 64;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 6;
    g.sigma[3] = 26;

    g.nc = 128;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 7;
    g.sigma[3] = 57;

    g.nc = 152;
  }

  if(g.coloring_distance == 8){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 36;
    g.sigma[3] = 52;

    g.nc = 256;
  }

  if(g.coloring_distance == 9){
    g.sigma[0] = 0;
    g.sigma[1] = 4;
    g.sigma[2] = 30;
    g.sigma[3] = 35;

    g.nc = 320;
  }

  if(g.coloring_distance == 10){
    g.sigma[0] = 0;
    g.sigma[1] = 8;
    g.sigma[2] = 31;
    g.sigma[3] = 60;

    g.nc = 481;
  }

  if(g.coloring_distance == 11){
    g.sigma[0] = 0;
    g.sigma[1] = 159;
    g.sigma[2] = 179;
    g.sigma[3] = 534;

    g.nc = 545;
  }

  if(g.coloring_distance == 12){
    g.sigma[0] = 0;
    g.sigma[1] = 10;
    g.sigma[2] = 54;
    g.sigma[3] = 55;

    g.nc = 704;
  }

  if(g.coloring_distance == 13){
    g.sigma[0] = 0;
    g.sigma[1] = 351;
    g.sigma[2] = 81;
    g.sigma[3] = 159;

    g.nc = 832;
  }

  if(g.coloring_distance == 14){
    g.sigma[0] = 0;
    g.sigma[1] = 443;
    g.sigma[2] = 363;
    g.sigma[3] = 325;

    g.nc = 1040;
  }

  if(g.coloring_distance == 15){
    g.sigma[0] = 0;
    g.sigma[1] = 585;
    g.sigma[2] = 359;
    g.sigma[3] = 759;

    g.nc = 1040;
  }

  if(g.coloring_distance == 16){
    g.sigma[0] = 0;
    g.sigma[1] = 1230;
    g.sigma[2] = 794;
    g.sigma[3] = 1730;

    g.nc = 1773;
  }

  if(g.coloring_distance == 17){
    g.sigma[0] = 0;
    g.sigma[1] = 64;
    g.sigma[2] = 413;
    g.sigma[3] = 65;

    g.nc = 2048;
  }

  if(g.coloring_distance == 18){
    g.sigma[0] = 0;
    g.sigma[1] = 2390;
    g.sigma[2] = 507;
    g.sigma[3] = 246;

    g.nc = 2624;
  }

  if(g.coloring_distance == 19){
    g.sigma[0] = 0;
    g.sigma[1] = 132;
    g.sigma[2] = 29;
    g.sigma[3] = 126;

    g.nc = 2774;
  }

  if(g.coloring_distance == 20){
    g.sigma[0] = 0;
    g.sigma[1] = 136;
    g.sigma[2] = 50;
    g.sigma[3] = 199;

    g.nc = 3200;
  }

  if(g.coloring_distance == 21){
    g.sigma[0] = 0;
    g.sigma[1] = 402;
    g.sigma[2] = 3905;
    g.sigma[3] = 64;

    g.nc = 4096;
  }

  if(g.coloring_distance == 22){
    g.sigma[0] = 0;
    g.sigma[1] = 64;
    g.sigma[2] = 364;
    g.sigma[3] = 1;

    g.nc = 4096;
  }

  if(g.coloring_distance == 23){
    g.sigma[0] = 0;
    g.sigma[1] = 3040;
    g.sigma[2] = 575;
    g.sigma[3] = 508;

    g.nc = 5120;
  }

  if(g.coloring_distance == 24){
    g.sigma[0] = 0;
    g.sigma[1] = 576;
    g.sigma[2] = 5506;
    g.sigma[3] = 1947;

    g.nc = 6144;
  }

  if(g.coloring_distance == 25){
    g.sigma[0] = 0;
    g.sigma[1] = 1237;
    g.sigma[2] = 477;
    g.sigma[3] = 1448;

    g.nc = 6356;
  }

  if(g.coloring_distance == 26){
    g.sigma[0] = 0;
    g.sigma[1] = 2144;
    g.sigma[2] = 2035;
    g.sigma[3] = 4048;

    g.nc = 7634;
  }

  if(g.coloring_distance == 27){
    g.sigma[0] = 0;
    g.sigma[1] = 92;
    g.sigma[2] = 2039;
    g.sigma[3] = 2564;

    g.nc = 8452;
  }

  if(g.coloring_distance == 28){
    g.sigma[0] = 0;
    g.sigma[1] = 271;
    g.sigma[2] = 418;
    g.sigma[3] = 30;

    g.nc = 8896;
  }

  if(g.coloring_distance == 29){
    g.sigma[0] = 0;
    g.sigma[1] = 455;
    g.sigma[2] = 134;
    g.sigma[3] = 21;

    g.nc = 9817;
  }

  if(g.coloring_distance == 30){
    g.sigma[0] = 0;
    g.sigma[1] = 4487;
    g.sigma[2] = 5327;
    g.sigma[3] = 1337;

    g.nc = 9941;
  }

  if(g.coloring_distance == 31){
    g.sigma[0] = 0;
    g.sigma[1] = 5219;
    g.sigma[2] = 1496;
    g.sigma[3] = 5805;

    g.nc = 13533;
  }

  if(g.coloring_distance == 32){
    g.sigma[0] = 0;
    g.sigma[1] = 5219;
    g.sigma[2] = 1496;
    g.sigma[3] = 5805;

    g.nc = 13533;
  }

  if(g.coloring_distance == 33){
    g.sigma[0] = 0;
    g.sigma[1] = 5219;
    g.sigma[2] = 1496;
    g.sigma[3] = 5805;

    g.nc = 13533;
  }

  if(g.coloring_distance == 34){
    g.sigma[0] = 0;
    g.sigma[1] = 269;
    g.sigma[2] = 302;
    g.sigma[3] = 608;

    g.nc = 19003;
  }

  if(g.coloring_distance == 35){
    g.sigma[0] = 0;
    g.sigma[1] = 835;
    g.sigma[2] = 4603;
    g.sigma[3] = 4704;

    g.nc = 22450;
  }

  if(g.coloring_distance == 36){
    g.sigma[0] = 0;
    g.sigma[1] = 835;
    g.sigma[2] = 4603;
    g.sigma[3] = 4704;

    g.nc = 22450;
  }

  if(g.coloring_distance == 37){
    g.sigma[0] = 0;
    g.sigma[1] = 835;
    g.sigma[2] = 4603;
    g.sigma[3] = 4704;

    g.nc = 22450;
  }

  if(g.coloring_distance == 38){
    g.sigma[0] = 0;
    g.sigma[1] = 2210;
    g.sigma[2] = 2244;
    g.sigma[3] = 98;

    g.nc = 25009;
  }

  if(g.coloring_distance == 39){
    g.sigma[0] = 0;
    g.sigma[1] = 2571;
    g.sigma[2] = 360;
    g.sigma[3] = 888;

    g.nc = 27632;
  }

  if(g.coloring_distance == 40){
    g.sigma[0] = 0;
    g.sigma[1] = 1154;
    g.sigma[2] = 1427;
    g.sigma[3] = 1023;

    g.nc = 31553;
  }

  if(g.coloring_distance == 41){
    g.sigma[0] = 0;
    g.sigma[1] = 6060;
    g.sigma[2] = 255;
    g.sigma[3] = 5234;

    g.nc = 33835;
  }

  if(g.coloring_distance == 42){
    g.sigma[0] = 0;
    g.sigma[1] = 3589;
    g.sigma[2] = 950;
    g.sigma[3] = 2125;

    g.nc = 34207;
  }

  if(g.coloring_distance == 43){
    g.sigma[0] = 0;
    g.sigma[1] = 1286;
    g.sigma[2] = 39;
    g.sigma[3] = 811;

    g.nc = 42466;
  }

  if(g.coloring_distance == 44){
    g.sigma[0] = 0;
    g.sigma[1] = 2615;
    g.sigma[2] = 659;
    g.sigma[3] = 2729;

    g.nc = 40993;
  }

  if(g.coloring_distance == 45){
    g.sigma[0] = 0;
    g.sigma[1] = 5321;
    g.sigma[2] = 5833;
    g.sigma[3] = 3641;

    g.nc = 44038;
  }

  if(g.coloring_distance == 46){
    g.sigma[0] = 0;
    g.sigma[1] = 5503;
    g.sigma[2] = 3367;
    g.sigma[3] = 1825;

    g.nc = 57332;
  }

  if(g.coloring_distance == 47){
    g.sigma[0] = 0;
    g.sigma[1] = 5503;
    g.sigma[2] = 3367;
    g.sigma[3] = 1825;

    g.nc = 57332;
  }

  if(g.coloring_distance == 48){
    g.sigma[0] = 0;
    g.sigma[1] = 1000;
    g.sigma[2] = 1759;
    g.sigma[3] = 68;

    g.nc = 61414;
  }

  if(g.coloring_distance == 49){
    g.sigma[0] = 0;
    g.sigma[1] = 3601;
    g.sigma[2] = 720;
    g.sigma[3] = 297;

    g.nc = 63283;
  }

  if(g.coloring_distance == 50){
    g.sigma[0] = 0;
    g.sigma[1] = 5857;
    g.sigma[2] = 1678;
    g.sigma[3] = 6073;

    g.nc = 70291;
  }

}

void dilution_check(){

  if(g.dilution != 1 && g.dilution != 2 && g.dilution != 3 && g.dilution != 4 && g.dilution != 12){
    printf("\nError: choose a correct dilution value (1, 2, 3, 4, 12)");
    exit(1);
  }

  if(g.dilution == 1)
    printf("\nNo dilution");

  if(g.dilution == 2)
    printf("\nPartial spin dilution");

  if(g.dilution == 3)
    printf("\nColor dilution");

  if(g.dilution == 4)
    printf("\nComplete spin dilution");

  if(g.dilution == 12)
    printf("\nSpin-Color dilution");
}

void graph_coloring(){

  MALLOC(g.num_colors, int, g.num_levels);
  MALLOC(g.dilution_ml, int, g.num_levels);

  if(g.my_rank==0){

    MALLOC(g.variances, double, g.num_levels);

    printf("\nProbing = %d\n", g.probing);
    printf("Coloring_distance = %d\n", g.coloring_distance);
    printf("Coloring_method = %d\n", g.coloring_method);
    
    if(g.probing_dimension == 3)
      get_sigma_3D();
    else
      get_sigma_4D();

    printf("sigma: %d %d %d %d\n", g.sigma[0], g.sigma[1], g.sigma[2], g.sigma[3]);
    printf("colors a the finest: %d\n", g.nc);

    dilution_check();

    double time_taken;

    double start_time = MPI_Wtime();

    g.colors = (int**)malloc(g.num_levels * sizeof(int*));

    if (g.colors == NULL)
        error0("Allocation error0\n");

    for(int level = 0; level < g.num_levels; level++){

      if(level == 0)
        g.dilution_ml[level] = g.dilution;
      else
        g.dilution_ml[level] = 1;

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

 /* 
    FILE *file = fopen("print_files/colors.txt", "w");

    for(int i = 0; i < g.num_levels; i++){
        fprintf(file, "\nColors at level %d\n [", i+1);

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
*/


  }
    MPI_Barrier(MPI_COMM_WORLD);
    setup_local_colors();
}
