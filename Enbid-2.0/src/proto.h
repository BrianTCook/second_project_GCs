void   allocate_memory(void);
void   begrun(void);
void   free_memory(void);
void   init(void);
void   read_ic0(char *fname);
void   read_ic1(char *fname);
void   read_ic2(char *fname);
void   read_parameter_file(char *fname);
void   savepositions(void );
void   savepositions_ioformat2(void);
void   savepositions_ioformat1(void);
void   savepositions_ioformat0(void);
void   set_sph_kernel(int TypeOfKernel, int dim);
void density_par(void);
size_t my_fwrite1(void *ptr,int* pattern,size_t nmemb,int flag_swap,FILE *stream);
size_t my_fread1(void *ptr, int *pattern, size_t nmemb,int flag_swap,FILE *stream);
void SwapEndian(void* addr, int* pattern);
template <class type> void order_particles(type* P,int *Id,int Nmax, int reverse);
void read_typelist(char *fname,vector<int> &npart, struct particle_data P[]);
void read_periodic_lengths(char *fname,double boxh[ND]);


