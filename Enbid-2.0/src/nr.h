/* Numerical recipies Shapiro, Teukolsky */
void nrerror(char error_text[]);
float *fvector(long nl, long nh);
void free_fvector(float *v, long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
void jacobi(float **a, int n, float d[], float **v, int *nrot);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);




