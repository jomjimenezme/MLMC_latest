#include "main.h"
//CURRENTLY WORKS ONLY FOR LATTICES WHERE T >= Z=Y=X
static const unsigned char RB_4D[16] = {
     0,  8,  9,  1,
    10,  2,  3, 11,
    12,  4,  5, 13,
     6, 14, 15,  7
};

static const unsigned char RB_3D[8] = {
    0, 4, 5, 1,
    6, 2, 3, 7
};

static const unsigned int perm_4D[16] = {
     0,  8,  4,  12,
    2,  10,  6, 14,
    1,  9,  5, 13,
     3, 11, 7,  15
};

uint32_t *g_hp_loc_4d = NULL;
uint32_t *g_hp_loc_3d = NULL;

// One-time precomputation of the hierarchical-probing bit string per local site.
// Bit p of g_hp_loc_*[site] == pi[p] of build_H / build_H_3d.
void hp_loc_setup( level_struct *l ){
  int *ll = l->local_lattice;
  int nsites = ll[0]*ll[1]*ll[2]*ll[3];
  int k0=g.global_k[0][0], k1=g.global_k[0][1];

  MALLOC( g_hp_loc_4d, uint32_t, nsites );
  MALLOC( g_hp_loc_3d, uint32_t, nsites );

  for( int s=0; s<nsites; s++ ){
    int lc[4], gc[4];
    lc[3] =  s % ll[3];
    lc[2] = (s /  ll[3]) % ll[2];
    lc[1] = (s / (ll[3]*ll[2])) % ll[1];
    lc[0] =  s / (ll[3]*ll[2]*ll[1]);
    for( int mu=0; mu<4; mu++ ) gc[mu] = g.my_coords[mu]*ll[mu] + lc[mu];

    // ---- 4-D packing: same bit logic as build_H ----
    uint32_t loc4 = 0; int c4 = 0;
    for( int k=0; k<k0; k++ ){
      if( k<k1 ){
        int dec = 8*((gc[0]>>k)&1) + 4*((gc[1]>>k)&1) + 2*((gc[2]>>k)&1) + ((gc[3]>>k)&1);
        dec = RB_4D[dec];
        for( int f=3; f>=0; f-- ){ loc4 |= (uint32_t)((dec>>f)&1) << c4; c4++; }
      } else {
        loc4 |= (uint32_t)((gc[0]>>k)&1) << c4; c4++;
      }
    }
    g_hp_loc_4d[s] = loc4;

    // ---- 3-D packing: same bit logic as build_H_3d ----
    uint32_t loc3 = 0; int c3 = 0;
    for( int k=0; k<k1; k++ ){
      int dec = 4*((gc[1]>>k)&1) + 2*((gc[2]>>k)&1) + ((gc[3]>>k)&1);
      dec = RB_3D[dec];
      for( int f=2; f>=0; f-- ){ loc3 |= (uint32_t)((dec>>f)&1) << c3; c3++; }
    }
    g_hp_loc_3d[s] = loc3;
  }

  #if 0  // one-time verification against reference implementation, remove after first successful run
  int *gl = g.global_lattice[0];
  for( int s=0; s<nsites; s++ ){
    int lc[4], gc[4];
    lc[3]=s%ll[3]; lc[2]=(s/ll[3])%ll[2]; lc[1]=(s/(ll[3]*ll[2]))%ll[1]; lc[0]=s/(ll[3]*ll[2]*ll[1]);
    for(int mu=0;mu<4;mu++) gc[mu]=g.my_coords[mu]*ll[mu]+lc[mu];
    int gidx = ((gc[0]*gl[1]+gc[1])*gl[2]+gc[2])*gl[3]+gc[3];
    for( int m=0; m<64; m++ ){
      int ref4 = build_H(gidx, m, 0);
      int new4 = (__builtin_popcount(g_hp_loc_4d[s] & (uint32_t)m)&1) ? -1 : +1;
      int ref3 = build_H_3d(gidx, m, 0);
      int new3 = (__builtin_popcount(g_hp_loc_3d[s] & (uint32_t)m)&1) ? -1 : +1;
      if( ref4!=new4 || ref3!=new3 )
        error0("hp_loc mismatch at site %d, m %d\n", s, m);
    }
  }
  printf0("hp_loc verification passed\n");
#endif

}

void hp_loc_free( level_struct *l ){
  int *ll = l->local_lattice;
  int nsites = ll[0]*ll[1]*ll[2]*ll[3];
  if(g_hp_loc_4d){ FREE( g_hp_loc_4d, uint32_t, nsites ); g_hp_loc_4d = NULL; }
  if(g_hp_loc_3d){ FREE( g_hp_loc_3d, uint32_t, nsites ); g_hp_loc_3d = NULL; }
}

//result[0]      = LSB  (least significative bit)
//result[bits-1] = MSB  (most significative bit)
int *dec2Bin(long long n, int bits){
  int *result = malloc(bits * sizeof(int));
  if (!result)
    return NULL;

  for (int i = 0; i < bits; i++){
    result[i] = n & 1;
    n >>= 1;
  }

  return result;
}

long long bin2Dec(int *bits_array, int bits) {
    long long result = 0;
    for (int i = 0; i < bits; i++) {
        result += bits_array[i] * (1LL << i);
    }
    return result;
}

void index_to_coord(int i, int coords[4], int level){
  coords[0] = i/(g.global_lattice[level][1]*g.global_lattice[level][2]*g.global_lattice[level][3]);
  coords[1] = (i/(g.global_lattice[level][2]*g.global_lattice[level][3]))%g.global_lattice[level][1];
  coords[2] = (i/g.global_lattice[level][3])%g.global_lattice[level][2];
  coords[3] = i%g.global_lattice[level][3];
}

void tensor_product_hadamard_dilution(const int *h, int N, const int *e, int *v){
  for(int site = 0; site < N; site++){
    for(int cs = 0; cs < 12; cs++){
      v[12*site + cs] = h[site] * e[cs];
    }
  }
}

int build_H(int i, int j, int level){

  int coords[4];
  index_to_coord(i, coords, level);

  int *t = dec2Bin(coords[0], g.global_k[level][0]);
  int *z = dec2Bin(coords[1], g.global_k[level][1]);
  int *y = dec2Bin(coords[2], g.global_k[level][2]);
  int *x = dec2Bin(coords[3], g.global_k[level][3]);

  int total_bits = g.global_k[level][0] + g.global_k[level][1] + g.global_k[level][2] + g.global_k[level][3];
  int *pi = malloc(total_bits * sizeof(int));

  int count = 0;
  for(int k=0; k<g.global_k[level][0]; k++){
    if(k<g.global_k[level][1]){
      int dec = 8*t[k] + 4*z[k] + 2*y[k] + x[k];
      dec = RB_4D[dec];
      int *rb = dec2Bin(dec, 4);
      for(int for_index = 3; for_index >= 0; for_index--) {
        pi[count] = rb[for_index];
        count++;
      }
      free(rb);
    }else{
       pi[count] = t[k];
       count++;
    }
  }

  int *pj = dec2Bin(j, total_bits);
/*
  if(g.anisotropic[level]==1){
    int pj_bits[4];
    int lsb = (g.k[level]-1)*4;

    pj_bits[0]=pj[lsb+1];
    pj_bits[1]=pj[lsb+2];
    pj_bits[2]=pj[lsb+3];
    pj_bits[3]=pj[lsb+4];

    int dec_pj_bits = 8*pj_bits[3] + 4*pj_bits[2] + 2*pj_bits[1] + pj_bits[0];
    int perm_j = perm_4D[dec_pj_bits];
    int *bin_perm_j = dec2Bin(perm_j, 4);

    pj[lsb+1] = bin_perm_j[0];
    pj[lsb+2] = bin_perm_j[1];
    pj[lsb+3] = bin_perm_j[2];
    pj[lsb+4] = bin_perm_j[3];

    free(bin_perm_j);

    int j_prime = bin2Dec(pj, total_bits);
    //if(g.my_rank==0) printf("\nOriginal j = %d -> permuted = %d", j, j_prime);
  }
*/
  int popcount = 0;
  for(int p=0; p<total_bits; p++)
    if(pi[p] == 1 && pj[p] == 1)
      popcount++;

  free(pi);
  free(pj);
  free(t);
  free(z);
  free(y);
  free(x);

  int h_ij;
  if(popcount%2==0)
    h_ij=1;
  else
    h_ij=-1;

  return h_ij;

}

int build_H_3d(int i, int j, int level){

  int coords[4];
  index_to_coord(i, coords, level);

  int *z = dec2Bin(coords[1], g.global_k[level][1]);
  int *y = dec2Bin(coords[2], g.global_k[level][2]);
  int *x = dec2Bin(coords[3], g.global_k[level][3]);

  int total_bits = g.global_k[level][1] + g.global_k[level][2] + g.global_k[level][3];
  int *pi = malloc(total_bits * sizeof(int));

  int count = 0;
  for(int k=0; k<g.global_k[level][1]; k++){
    int dec = 4*z[k] + 2*y[k] + x[k];
    dec = RB_3D[dec];
    int *rb = dec2Bin(dec, 3);
    for(int for_index = 2; for_index >= 0; for_index--){
      pi[count] = rb[for_index];
      count++;
    }
    free(rb);
  }

  int *pj = dec2Bin(j, total_bits);
  int popcount = 0;
  for(int p=0; p<total_bits; p++)
    if(pi[p] == 1 && pj[p] == 1)
      popcount++;

  free(pi);
  free(pj);
  free(z);
  free(y);
  free(x);

  int h_ij;
  if(popcount%2==0)
    h_ij=1;
  else
    h_ij=-1;

  return h_ij;

}


