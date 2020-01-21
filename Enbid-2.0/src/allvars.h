/* This file declares all global variables. Further variables should be added here,
   and declared as 'extern'. The actual existence of these variables is provided by
   the file 'allvars.cpp'. To produce 'allvars.cpp' from 'allvars.h', do the following:

   1.) Erase all #define's
   2.) add #include "allvars.h"
   3.) delete all keywords 'extern'
   4.) delete all struct definitions enclosed in {...}, e.g.
   "extern struct global_data_all_processes {....} All;"
   becomes "struct global_data_all_processes All;"
   5.) delete extern struct xyz_data {...}; completely
*/

#include <stdio.h>
#include<iostream>
#include<vector>
using namespace std;

#ifdef MEDIAN
#define ROOT		1
#define LOWER(i,node)	(i<<1)
#define UPPER(i,node)	((i<<1)+1)
#define PARENT(i,node)	(i>>1)
#define SIBLING(i,node) 	((i&1)?i-1:i+1)
#define SETNEXT(i,node)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}
#else
#define ROOT	1
#define LOWER(i,node)	(node[i].left)
#define UPPER(i,node)	(node[i].left+1)
#define PARENT(i,node)	(node[i].parent)
#define SIBLING(i,node) 	((i&1)?i-1:i+1)
#define SETNEXT(i,node)\
{\
	while (i&1) i=node[i].parent;\
	++i;\
	}
#endif


#ifdef DIM6
#define ND 6  //const int ND=6;
#endif
#ifdef DIM3
#define ND 3  //const int ND=3;
#endif
#ifdef DIMO
#define ND 4 //const int ND=5;
#endif

#ifdef T3E
typedef short int int4byte;   /* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
#endif

#define real double
#define real1 float

/* ... often used constants (cgs units) */
#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37
#define  THIRD            (1.0/3.0)
#ifndef  PI
#define  PI               3.14159265358979323846
#endif
#define  PI_INV           (1/PI)
#define  LN2              0.69314718
#define  KERNEL_TABLE 1000
#define  MAX_NGB  1000000  /* 20000defines maximum length of neighbour list */

//extern int    BatchFlag;
extern int    NumPart;

/* variables for input/output ,  usually only used on process 0 */
extern  char   ParameterFile[100];

/* tabulated smoothing kernel */
extern double  Kernel[KERNEL_TABLE+2],
    KernelDer[KERNEL_TABLE+2],
    KernelDer2[KERNEL_TABLE+2],
    KernelRad[KERNEL_TABLE+2];

/* this struct contains mostly code parameters read from the parameter file */
extern struct global_data_all_processes
{
    /* Code options */
    int    TypeOfSmoothing;
    int    KernelBiasCorrection;
    int    AnisotropicKernel;
    int    Dimensions;
    int    VolCorr;
    int    TypeOfKernel;
    double SpatialScale;
    int    PartBoundary;
    int    NodeSplittingCriterion;
    int    CubicCells;
    double Anisotropy;
    int    DesNumNgb;
    int    NumBucket;
    int    DesNumNgbA;
    int    NumBucketA;
    /* Cosmology */
//    double BoxSize, BoxHalf;
    int    MedianSplittingOn;
    float   hs[ND],hsv;
    double  MassTable[6];
    /* File options */
    int   ICFormat;
    double OmegaBaryon;
    double Omega0;
    char    InitCondFile[100],SnapshotFileBase[100];
    char   InputDir[100],
	InitFileBase[100];
    /* Some other global parameters */
    int   MaxPart,TypeListOn;
    int flag_swap;
    int order_flag;
    int PeriodicBoundaryOn;
#ifdef PERIODIC
    double boxh[ND];
#endif

} All;




extern struct particle_data
{
    int  ID,Type;//,NumNgb;           /* unique particle identifier */
    float     Pos[ND];       // particle position
#ifdef DIM3
    float     Vel[3];       // particle velocity
#endif
    float     Mass;         /* particle mass */
    float     Density;
} *P,*P_data,*Part;//,*P1;


extern vector<int> npart,npartc; 

/* Header for the standard file format. */
extern struct io_header_1
{
    int4byte npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int4byte flag_sfr;
    int4byte flag_feedback;
    int4byte npartTotal[6];
    int4byte flag_cooling;
    int4byte num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    /* extra flags other than gadget*/
    int4byte flag_id;          /* if IDs needed in output */
    int4byte flag_dim;         /* no of dimensions */
    /* signifies that output file contains estimated density */
    int4byte flag_density;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8-3*4];  /* fills to 256 Bytes */
} header1;

extern double t00,t11,t22;
extern double xbmax[ND],xbmin[ND],kbuf[ND]; /* boundary calculation variables */


/* extern struct param_data */
/* { */
/*     /\* Code options *\/ */
/*     int    BeginFileNum,EndFileNum,FileNum; */
/* } Var; */





extern int *list_kdt,*list_kd;//,*father;
extern float *listr;
extern real searchx[ND][2], searchcenter[ND],metric[ND],metric1[ND],gmatrix[ND][ND];//,h[ND];
extern vector<int> trees;
extern FILE *treefile_fpt;
//extern real BoxHalf, Box;
// count = no of particles
// k1=dim to be split  ,kcut= split location
// for leaf node k1 is index to its pnode
// left and right left= index(pointer) to left and right nodes

#ifdef MEDIAN
extern struct NODE
{
    int count,k1,lid;  // 2*20=40 bytes,
    real kcut;
    struct XNODE * bnd;
} *nodes,*nodes_base,*nodesC;
#else
extern struct NODE
{
    real kcut;
    int count,k1,lid,left,parent;  // 2*20=40 bytes,
    struct XNODE * bnd;
} *nodes,*nodes_base,*nodesC;
#endif

//contains additional details of node
// lid: first id of particle
// k0:  last dimension split
// bleft: 0 or 1 checks if the surface is derived from boundary or not
extern struct SNODE
{
    //  real Density,h[ND];//,ex,ev,e[ND];  2*20=40bytes
    int  lid,k0;
    bool bleft[ND],bright[ND];
} *snodes,*snodes_base;

// contains co-ordinates (surfaces) of leaf node
extern struct PNODE
{
    real x[ND][2];              // 48 bytes
} *pnodes,*pnodes_base;


// contains co-ordinates (surfaces) of all nodes
extern struct XNODE
{
    real x[ND][2];            // 2*48=96 bytes +12bytes =
} *xnodes,*xnodes_base;



extern struct linklist 
{
    real r;
    int p;
    float x[ND];
} *pqHead;

extern struct pqueue 
{
    struct linklist* pq;
    bool operator<( const struct pqueue &a ) const
    {
	return pq->r < a.pq->r;
    }
    bool operator>( const struct pqueue &a ) const
    {
	return pq->r > a.pq->r;
    }
}*pqx1;


extern float **mrho,*d,**ve;
extern int pnext,pnew;




