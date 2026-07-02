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
      for(int for_index=0; for_index<4; for_index++){
        pi[count]=rb[for_index];
	count++;
      }
      free(rb);
    }else{
       pi[count] = t[k];
       count++;
    }
  }

  int *pj = dec2Bin(j, total_bits);

  if(g.anisotropic[level]==1){
    int pj_bits[4];
    int lsb = (g.k[level]-2)*4;

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
    if(g.my_rank==0) printf("\nOriginal j = %d -> permuted = %d", j, j_prime);
  }

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
    for(int for_index=0; for_index<3; for_index++){
      pi[count]=rb[for_index];
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


