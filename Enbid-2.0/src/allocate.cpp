#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "functions.h"


/* This routine allocates memory for 
 * particle storage.
 */
void allocate_memory(void)
{
  
  
    int bytes=0,bytes_tot=0;
  
    if(All.MaxPart>0)
    {
	if(!(P_data=new particle_data[All.MaxPart]))
	{
	    printf("failed to allocate memory for `P_data' (%d bytes).\n",bytes);
	    endrun(1);
	}
	bytes_tot+=All.MaxPart*sizeof(struct particle_data);
    
	P= P_data-1;   /* start with offset 1 */
	Part= &P[1];
	printf("Allocated %g MByte for particle storage.\n",bytes_tot/(1024.0*1024.0));
    }
  
}



/* This routine frees the memory for the particle storage,
 */
void free_memory(void)
{
  
    if(All.MaxPart>0)
    {
	free(P_data);
    }
  
}
