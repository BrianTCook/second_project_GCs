#include "allvars.h"
//int    BatchFlag;
int    NumPart;

/* variables for input/output ,  usually only used on process 0 */
char   ParameterFile[100];

/* tabulated smoothing kernel */
double  Kernel[KERNEL_TABLE+2],
    KernelDer[KERNEL_TABLE+2],
    KernelDer2[KERNEL_TABLE+2],
    KernelRad[KERNEL_TABLE+2];

/* this struct contains mostly code parameters read from the parameter file */
struct global_data_all_processes All;

vector<int> trees;
FILE *treefile_fpt;
//real BoxHalf, Box;

/* The following structure holds all the information that is
 * stored for each particle.
 */

vector<int> npart,npartc;

struct io_header_1 header1;

//struct param_data Var;

struct particle_data  *P,*P_data,*Part;


/* Header for the GADGET file format. */

double t00,t11,t22;
double xbmax[ND],xbmin[ND],kbuf[ND]; /* boundary calculation variables */


real searchx[ND][2], searchcenter[ND],metric[ND],metric1[ND],gmatrix[ND][ND];
int *list_kdt,*list_kd;
float *listr;

struct NODE
*nodes,*nodes_base,*nodesC;
struct SNODE
*snodes,*snodes_base;
struct PNODE
*pnodes,*pnodes_base;


struct XNODE
*xnodes,*xnodes_base;

struct linklist  *pqHead;
struct pqueue *pqx1;
int pnext,pnew;

float **mrho,*d,**ve;



