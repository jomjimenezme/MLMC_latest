#include "main.h"
#include "data_layout.h"
#include "stdbool.h"

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

int min(int v[], int size) {
    int min = v[0];

    for (int i = 1; i < size; i++) {
        if (v[i] < min) {
            min = v[i];
        }
    }

    return min;
}

int pow_int(int base, int exp) {
    int res = 1;
    for (int i = 0; i < exp; ++i) {
        res *= base;
    }
    return res;
}

int pow2_valuation(unsigned int x) {
    return __builtin_ctz(x);
}

void vector_copy(int *dest, int *src, int size) {
    for (int i = 0; i < size; i++) {
        dest[i] = src[i];
    }
}

bool contains(int *array, int size, int value) {
    for (int i = 0; i < size; i++) {
        if (array[i] == value)
            return true;
    }
    return false;
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


void get_coloring_dimension(){

  if(g.trace_op_type == 1 || g.trace_op_type == 2 || g.trace_op_type == 6 || g.trace_op_type == 7 || g.trace_op_type == 8)
    g.probing_dimension = 3;

  if(g.trace_op_type == 9 || g.trace_op_type == 10 || g.trace_op_type == 11 || g.trace_op_type == 12)
    g.probing_dimension = 4;
}

void setup_local_colors_direct(){
    int dims[4], periods[4], coords[4];
    MPI_Cart_get(g.comm_cart, 4, dims, periods, coords);

    g.local_colors = (int**)malloc(g.num_levels * sizeof(int*));

    for(int level = 0; level < g.num_levels; level++){
        int T = g.global_lattice[level][0];
        int Z = g.global_lattice[level][1];
        int Y = g.global_lattice[level][2];
        int X = g.global_lattice[level][3];
        int Nt_loc = T / dims[0];
        int Nz_loc = Z / dims[1];
        int Ny_loc = Y / dims[2];
        int Nx_loc = X / dims[3];
        const int local_size = Nt_loc * Nz_loc * Ny_loc * Nx_loc;

        MALLOC(g.local_colors[level], int, local_size);

        // global coordinate offsets of this rank's block
        int t0 = coords[0] * Nt_loc;
        int z0 = coords[1] * Nz_loc;
        int y0 = coords[2] * Ny_loc;
        int x0 = coords[3] * Nx_loc;

        int idx = 0;
        for (int lt = 0; lt < Nt_loc; lt++)
            for (int lz = 0; lz < Nz_loc; lz++)
                for (int ly = 0; ly < Ny_loc; ly++)
                    for (int lx = 0; lx < Nx_loc; lx++, idx++)
                    {
                        if( level == 0 ){
                            int t = t0 + lt;
                            int z = z0 + lz;
                            int y = y0 + ly;
                            int x = x0 + lx;
                            int col = ( g.sigma[0]*t + g.sigma[1]*z
                                      + g.sigma[2]*y + g.sigma[3]*x ) % g.num_colors[0];
                            // manual two-form override for d=7 would replace the line above:
                            //   96^3: int col = 3*(( z + 2*y + 11*x ) % 48) + ( y % 3 );
                            //   64^3: int col = 2*(( z + 3*y + 13*x ) % 64) + (( y + x ) % 2);
                            g.local_colors[level][idx] = col + 1;   // 1-based
                        } else {
                            g.local_colors[level][idx] = 1;         // coarse levels: single color
                        }
                    }
    }
}

void get_sigma_4D(){
  
  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;
    
    g.num_colors[0] = 2;
  }
  
  if(g.coloring_distance == 2){
    g.sigma[0] = 1;
    g.sigma[1] = 2;
    g.sigma[2] = 3;
    g.sigma[3] = 4;
    
    g.num_colors[0] = 10;
  }
  
  if(g.coloring_distance == 3){
    g.sigma[0] = 1;
    g.sigma[1] = 5;
    g.sigma[2] = 55;
    g.sigma[3] = 61;
    
    g.num_colors[0] = 16;
  }
  
  if(g.coloring_distance == 4){
    g.sigma[0] = 1;
    g.sigma[1] = 8;
    g.sigma[2] = 12;
    g.sigma[3] = 18;
    
    g.num_colors[0] = 64;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 38;
    g.sigma[1] = 1;
    g.sigma[2] = 12;
    g.sigma[3] = 16;

    g.num_colors[0] = 128;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 3;
    g.sigma[1] = 20;
    g.sigma[2] = 48;
    g.sigma[3] = 50;

    g.num_colors[0] = 320;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 40;
    g.sigma[1] = 32;
    g.sigma[2] = 33;
    g.sigma[3] = 61;

    g.num_colors[0] = 512;
  }
  
}

void get_sigma_3D(){
  
  if(g.coloring_distance == 1){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;
    
    g.num_colors[0] = 2;
  }
  
  if(g.coloring_distance == 2){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 2;
    g.sigma[3] = 3;
    
    g.num_colors[0] = 8;
  }
  
  if(g.coloring_distance == 3){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 3;
    g.sigma[3] = 5;
    
    g.num_colors[0] = 16;
  }
  
  if(g.coloring_distance == 4){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 6;
    g.sigma[3] = 9;
    
    g.num_colors[0] = 32;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 11;
    g.sigma[3] = 27;

    g.num_colors[0] = 88;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 8;
    g.sigma[3] = 44;

    g.num_colors[0] = 128;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 9;
    g.sigma[3] = 33;

    g.num_colors[0] = 176;
  }

  if(g.coloring_distance == 8){
    g.sigma[0] = 0;
    g.sigma[1] = 7;
    g.sigma[2] = 48;
    g.sigma[3] = 51;

    g.num_colors[0] = 272;
  }

  if(g.coloring_distance == 9){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 33;
    g.sigma[3] = 45;

    g.num_colors[0] = 352;
  }
  
}

void dilution_check(int level){

  if(g.dilution[level] != 1 && g.dilution[level] != 2 && g.dilution[level] != 3 && g.dilution[level] != 4 && g.dilution[level] != 12){
    printf("\nError: choose a correct dilution value (1, 2, 3, 4, 12) at level %d", level);
    exit(1);
  }

  if(g.dilution[level] == 1)
    if(g.my_rank==0) printf("\nNo dilution at level %d\n", level);

  if(g.dilution[level] == 2)
    if(g.my_rank==0) printf("\nPartial spin dilution at level %d\n", level);

  if(g.dilution[level] == 3)
    if(g.my_rank==0) printf("\nColor dilution at level %d\n", level);

  if(g.dilution[level] == 4)
    if(g.my_rank==0) printf("\nComplete spin dilution at level %d\n", level);

  if(g.dilution[level] == 12)
    if(g.my_rank==0) printf("\nSpin-Color dilution at level %d\n", level);
}

void coloring_scheme(){
  // sigma / num_colors selection: pure table lookup, must now run on ALL ranks
  // (previously rank-0 only; other ranks learned the coloring via the Bcast)
  if(g.probing_dimension == 3)
    get_sigma_3D();
  else
    get_sigma_4D();

  // legacy safety Bcasts (harmless if already consistent on all ranks)
  MPI_Bcast(g.num_colors, g.num_levels, MPI_INT, 0, g.comm_cart);
  MPI_Bcast(g.dilution,   g.num_levels, MPI_INT, 0, g.comm_cart);

  if(g.my_rank==0){
    printf("\nProbing = %d - Classical probing\n", g.probing);
    printf("Coloring_distance = %d\n", g.coloring_distance);
    printf("Grids to be colored = %d\n", g.colored_grids);
    printf("Coloring dimension = %d\n", g.probing_dimension);
    printf("sigma: %d %d %d %d\n", g.sigma[0], g.sigma[1], g.sigma[2], g.sigma[3]);
  }

  // torus validity: sigma_mu * n_mu must vanish mod nc, otherwise the coloring
  // silently violates its distance across the periodic boundary
  for( int mu=0; mu<4; mu++ )
    if( ( g.sigma[mu] * g.global_lattice[0][mu] ) % g.num_colors[0] != 0 )
      error0("multiplier coloring INVALID on this torus: sigma[%d]=%d, n=%d, nc=%d\n",
             mu, g.sigma[mu], g.global_lattice[0][mu], g.num_colors[0]);

  for( int level = 1; level < g.num_levels; level++ )
    g.num_colors[level] = 1;

  for( int level = 0; level < g.num_levels; level++ )
    dilution_check(level);

  setup_local_colors_direct();

  if(g.my_rank==0)
    for( int level = 0; level < g.num_levels; level++ )
      printf("\n Colors at depth %d : \t %d \n", level, g.num_colors[level]);

  MPI_Barrier(MPI_COMM_WORLD);
}

int* find_indices(int *array, int size, int value, int *count){
    // first pass: count how many matches
    *count = 0;
    for(int i = 0; i < size; i++){
        if(array[i] == value)
            (*count)++;
    }
    // allocate result
    int *indices = (int*)malloc(*count * sizeof(int));
    if(indices == NULL)
        return NULL;
    // second pass: fill indices
    int j = 0;
    for(int i = 0; i < size; i++){
        if(array[i] == value)
            indices[j++] = i;
    }
    return indices;
}

void graph_coloring(){

  get_coloring_dimension();
  if(g.probing == 1) coloring_scheme();

 if(g.probing == 2){
   if(g.my_rank==0) printf("\nProbing = %d - Hierarchical probing\n", g.probing);
   if(g.my_rank==0) printf("Applied to num levels = %d\n", g.colored_grids);
   if(g.my_rank==0) printf("Probing dimension = %d\n", g.probing_dimension);
   for(int level = 0; level<g.num_levels; level++){
     g.global_k[level][0] = pow2_valuation(g.global_lattice[level][0]);
     g.global_k[level][1] = pow2_valuation(g.global_lattice[level][1]);
     g.global_k[level][2] = pow2_valuation(g.global_lattice[level][2]);
     g.global_k[level][3] = pow2_valuation(g.global_lattice[level][3]);
     dilution_check(level);
     g.num_colors[level] = g.n_had[level];

     if(g.my_rank==0){
       printf("Number of Hadamard vectors at level %d = %d\n", level, g.num_colors[level]);
       printf("Global k_t at level %d = %d\n", level, g.global_k[level][0]);
       printf("Global k_z at level %d = %d\n", level, g.global_k[level][1]);
       printf("Global k_y at level %d = %d\n", level, g.global_k[level][2]);
       printf("Global k_x at level %d = %d\n", level, g.global_k[level][3]);
       printf("Hierarchy at level %d = %d\n", level, g.k[level]);
       printf("Anisotropy at level %d = %d\n", level, g.anisotropic[level]);
     }
   }
//MOVE THIS CHECK INSIDE THE PER-LEVEL LOOP WHEN DOING MGMLMC
   int tb = (g.probing_dimension == 4)
             ? g.global_k[0][0] + 3*g.global_k[0][1]
             : 3*g.global_k[0][1];
   long long limit = 3LL * (1LL << tb);
   if( (long long)g.num_colors[0] > limit )
     error0("HP: num_colors %d exceeds 3*2^%d = %lld (one three-coloring level past the bit levels); undefined beyond\n",
             g.num_colors[0], tb, limit);

 }
 
 if(g.probing == 0){
   if(g.my_rank==0) printf("\nProbing = %d - Plain estimator\n", g.probing);
      for(int level = 0; level<g.num_levels; level++){
         g.global_k[level][0] = pow2_valuation(g.global_lattice[level][0]);
         g.global_k[level][1] = pow2_valuation(g.global_lattice[level][1]);
         g.global_k[level][2] = pow2_valuation(g.global_lattice[level][2]);
         g.global_k[level][3] = pow2_valuation(g.global_lattice[level][3]);
         dilution_check(level);
         g.num_colors[level] = 1;
   }
   g.probing = 2;
 }
}
