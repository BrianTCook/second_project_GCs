#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "tree.h"
#include "functions.h"

/*
 *  init() reads in the initial conditions,
 *  and allocates storage for the tree(s).
 *  Then generates the tree  and finally calculates
 *  the density.
 */
void init(void)
{
    int i,j;
    float h2,m_min,m_max;
    double t1,t2,xmean[ND];

    t11=second();
//     if(BatchFlag==1)
// 	sprintf(All.InitCondFile,"%s%s_%03d",All.InputDir,All.InitFileBase,Var.FileNum);

    if(All.AnisotropicKernel==1)
    {
	if((All.TypeOfSmoothing!=2)and(All.TypeOfSmoothing!=3))
	{
	    cout<<"For AnisotropicKernel=1, TypeOfSmoothing should be 2 or 3"<<endl; 
	    endrun(10);
	}
    }


    switch(All.ICFormat)
    {
	case 0:read_ic0(All.InitCondFile);break;
	case 1:read_ic1(All.InitCondFile);break;
	case 2:read_ic2(All.InitCondFile);break;
	default:cout<<"Parameter ICFormat should be between 0 to 2"<<endl;endrun(10);
    }


    if(All.TypeListOn)
    {
	char buf[100];
	char file_suf[10]="_typelist";
	sprintf(buf,"%s%s",All.InitCondFile,file_suf);
	read_typelist(buf,npart,P);	    
    }

#ifdef PERIODIC
    if(All.PeriodicBoundaryOn==1)
    {
	char file_periodic[21]="periodic_lengths.txt";
	read_periodic_lengths(file_periodic,All.boxh);	    

	for(j=0,i=0;j<int(npart.size());j++)
	    if(npart[j]>0) i++;
	if(i>1)
	{
	    cout<<"WARNING:multiple species present,"<<endl; 
#ifdef WARN
	    checkrun("a constant global periodic boundary will be used for all of them","exit","continue");
#else
	    cout<<"a constant global periodic boundary will be used for all of them"<<endl;
#endif
	}

    }
    else
    {
	cout<<"Code compiled with PERIODIC option, Set parameter PeriodicBoundaryOn to 1"<<endl;
	endrun(10);
    }
#else
    if(All.PeriodicBoundaryOn!=0)
    {
	cout<<"Set parameter PeriodicBoundaryOn to 0"<<endl;
	endrun(10);
    }
#endif

// check multiple mass and issue warning
    for(j=0;j<int(npart.size());j++)
    {
	i=npartc[j]+1;
	m_min=P[i].Mass; m_max=P[i].Mass;
	for(i=npartc[j]+1;i<=(npartc[j]+npart[j]);i++)
	    {
		if(m_min>P[i].Mass) m_min=P[i].Mass;
		if(m_max<P[i].Mass) m_max=P[i].Mass;
	    }
	if(m_min!=m_max)
	{
	    if(All.TypeOfSmoothing==1)
	    {
		cout<<"For FiEstAS smoothing multiple mass support for same species of particle is not available"<<endl;
		endrun(10);
		}
	    else
	    {
#ifdef WARN
		checkrun("WARNING! particles with multiple mass results might not be accurate","exit","continue");
#else
		cout<<"WARNING! particles with multiple mass results might not be accurate"<<endl;
#endif
	    }
	    
	}
    }
    


    header1.flag_density=1;
    treeallocate( All.MaxPart);


    /* scaling the velocity and spatial co-ordinates */
    for(j=0;j<ND;j++)
	All.hs[j]=1.0;

    if(All.SpatialScale>0)
    {
	if(ND==6)
	{
	    for(j=0;j<3;j++)
		All.hs[j]=All.SpatialScale;
	}
    }
    else
    {
	for(j=0;j<ND;j++)
	    xmean[j]=0.0;
	for(i=1; i<=NumPart; i++) /*  start-up initialization */
	{
	    for(j=0;j<ND;j++)
	    {
		xmean[j]+=P[i].Pos[j];
	    }
	}
	for(j=0;j<ND;j++)
	    xmean[j]/=NumPart;

	for(i=1; i<=NumPart; i++) /*  start-up initialization */
	{
	    for(j=0;j<ND;j++)
	    {
		All.hs[j]+=(P[i].Pos[j]-xmean[j])*(P[i].Pos[j]-xmean[j]);
	    }
	}

	for(j=0;j<ND;j++)
	    All.hs[j]=sqrt(All.hs[j]/NumPart);
	cout<<"Scaling Co-ordinates based on variance"<<endl;
	h2=a_max(All.hs,0,ND-1);
	for(j=0;j<ND;j++)
	{
	    cout<<"i = "<<j<<" Sigma = "<<All.hs[j]<<" Scale = "<<All.hs[j]/h2<<endl;
	    All.hs[j]=All.hs[j]/h2;
	}

    }

    for(j=0,All.hsv=1.0;j<ND;j++)
	All.hsv*=All.hs[j];
    cout<<"Scaling Co-ordinates as x[i]=x[i]/h[i] with h[i]->"<<endl;
    for(j=0;j<ND;j++)
    {
	cout<<All.hs[j]<<" ";
    }
    cout<<"\n"<<endl;



    for(i=1; i<=NumPart; i++) /*  start-up initialization */
    {
	for(j=0;j<ND;j++)
	{
	    P[i].Pos[j]=P[i].Pos[j]/All.hs[j];
	}
    }

#ifdef PERIODIC
	for(j=0;j<ND;j++)
	{
	    if(All.boxh[j]>0)
	    {
		All.boxh[j]=All.boxh[j]/All.hs[j];
	    }
	}
#endif



    /* Build Tree */
    if((All.NodeSplittingCriterion==1)&&(All.CubicCells==1))
    {
	if((ND!=6)&&(ND!=3))
	{
	    cout<<"Cubic Cells not allowed for dimensions other than 3 or 6 here Dim="<<ND<<endl;
	    cout<<"Choose another option or modify routine"<<endl;
	    endrun(10);
	}
    }

    printf("Starting to Build Tree ......."); fflush(stdout);
    for(i=1; i<=NumPart; i++)
    {
	P[i].Density=0.0;
    }




    t1=second();
    treebuild_all();
    t2=second();
    cout<<"Treebuild time = "<<timediff(t1,t2)<<" s \n"<<endl;
    t22=second();  t00+=timediff(t11,t22);


    All.order_flag=1;
    if(All.order_flag==1)
    {
	for(i=0; i<NumPart; i++)
	    list_kdt[i]=list_kd[i];
	order_particles(P+1,list_kd,NumPart,0);
	for(i=0; i<NumPart; i++)
	    list_kd[i]=i;
    }




	density_par();


    if(All.order_flag==1)
	order_particles(P+1,list_kdt,NumPart,1);

//     for(int ii=0;ii<100;ii++)
// 	cout<<"problem "<<Part[ii].Density<<endl;

	savepositions();


    /* correct for global scaling */
    t22=second();  t00+=timediff(t11,t22);

}



template <class type> void order_particles(type* P,int *Id,int Nmax, int reverse)
{
    int i;
    type Psave;
    int idsource, idest,*Id1;

    Id1= (int *) malloc(sizeof(int)*Nmax);

    if(reverse==1)
	for(i = 0; i < Nmax; i++)
	    Id1[Id[i]]=i;
    else
	for(i = 0; i < Nmax; i++)
	    Id1[i]=Id[i];

    for(i = 0; i < Nmax; i++)
    {
	if(Id1[i] != i)
        {
            Psave=P[i];
            idest = i;
	    do
            {
		idsource=Id1[idest];
		if(idsource == i)
		{
		    P[idest]=Psave;
		    Id1[idest] =idest;
		    break;
		}
		P[idest] = P[idsource];
		Id1[idest] = idest;
		idest = idsource;
            }
	    while(1);
        }
    }
    free(Id1);
}





