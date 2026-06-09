int *dec2Bin(long long n, int bits);
long long bin2Dec(int *bits_array, int bits);
void index_to_coord(int i, int coords[4], int level);
int build_H(int i, int j, int level);
void tensor_product_hadamard_dilution(const int *h, int N, const int *e, int *v);
