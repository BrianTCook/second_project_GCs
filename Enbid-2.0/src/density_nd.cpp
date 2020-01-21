#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <algorithm>
#include "allvars.h"
#include "nr.h"
#include "tree.h"
#include "functions.h"
#include "ngb_search.h"
#include "proto.h"

void density_par(void)
{
    int i,j,ii,ji,pi=1,i1;
    double t0;
    struct NODE* nodes1=NULL;

    struct linklist *pqStart=NULL;
    struct pqueue *pqx=NULL;
    bool *imark=NULL,*idone=NULL;
    struct linklist *pqStartA=NULL;
    struct pqueue *pqxA=NULL;
    bool *imarkA=NULL,*idoneA=NULL;


    create_linklist(pqx,pqStart,imark,idone,All.DesNumNgb,NumPart);
    if(All.AnisotropicKernel==1)
    {
	d=fvector(1,ND);  ve=matrix(1,ND,1,ND);  mrho=matrix(1,ND,1,ND);
	create_linklist(pqxA,pqStartA,imarkA,idoneA,All.DesNumNgbA,NumPart);
    }



    cout<<"Density Calculation. Smoothing ....."<<endl;
    t0=t00;
    t22=second();


    t11=t22; t22=second();t00+=timediff(t11,t22);
    t0=t00;
    t22=second();
    /* Main loop over all particles */
    i=0;i1=1;pi=NumPart/100;

    for(j=0;j<ND;j++)
	metric[j]=1.0;


    for(j=0; j<int(npart.size()); j++)
    {

	i=npartc[j];
	if(npart[j]>1)
	{


//	    cout<<"Smoothing Particle Type "<<j<<" "<<npart[j]<<endl;
	    if(npart[j]<All.DesNumNgb) 
	    {
		cout<<"Min "<<All.DesNumNgb<<" particles needed for smoothing"<<endl;
		endrun(10);
	    }

	    nodes1=nodes+trees[j]-1;
	    if(All.AnisotropicKernel==1)
		initialize_linklist(pqStartA,imarkA,i,All.DesNumNgbA,NumPart);
	    initialize_linklist(pqStart,imark,i,All.DesNumNgb,NumPart);


	    pnew=1;
	    pnext=i+1;

	    for(ji=0; ji<npart[j]; ji++)
	    {
		ii=i;

 		while(idone[ii]==0)
 		{
 		    idone[ii]=1;
  		    Part[ii].Density= density_general(Part[ii].Pos, nodes1,All.DesNumNgbA,pqxA,pqStartA,idoneA,imarkA,All.DesNumNgb,pqx,pqStart,idone,imark);

		    if(i1==pi)
		    {
			pi=pi+NumPart/100;
			t11=t22; t22=second();t00+=timediff(t11,t22);
			printf("Evaluated = %3d %c Time Left = %f s of %f s Par no = %d Density = %e \n",((i1)*100)/NumPart,'%',(t00-t0)*(NumPart-i1-1)*1.0/(i1),(t00-t0)*NumPart*1.0/(i1),list_kdt[ii]+1,Part[ii].Density);		
			fflush(stdout);
		    }
		    i1++;

		    ii=pnext;
 		}
		i++;
		
	    }

	}
    }


//   	  order_particles(P+1,list_kdt,NumPart,1);
//   		for(i=0; i<NumPart; i=i+1000)
//   		    cout<<i<<" "<<Part[i].Density<<endl;

    if(All.AnisotropicKernel==1)
    {
	free_fvector(d, 1,ND);
	free_matrix(ve,1,ND,1,ND);
	free_matrix(mrho,1,ND,1,ND);
	delete [] pqStartA;
	delete [] pqxA;
	delete [] imarkA;
	delete [] idoneA;
    }
    delete [] pqStart;
    delete [] pqx;
    delete [] imark;
    delete [] idone;

    t11=t22; t22=second();t00+=timediff(t11,t22);
    cout<<"\nTotal Smoothing Time = "<<t00-t0<<" s"<<endl;
    cout<<All.VolCorr<<" "<<All.TypeOfSmoothing<<endl;
}





