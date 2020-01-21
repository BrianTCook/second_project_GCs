/***************************************************************************
                          tree.cpp  -  description
                             -------------------
    begin                : Mon Jan 16 2006
    copyright            : (C) 2006 by Sanjib Sharma
    email                : ssharma@aip.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "tree.h"
#include "functions.h"

#define NULL_N   -1
// #define real double
// #define real1 float


static int numnodestree[150],NumPart1;    /* number of (internal) nodes in each tree */
static int numnodestotal,nfree1;      /* total number of internal nodes */
static int MaxNodes;           /* maximum allowed number of internal nodes */
static int nnodes,pi,rootnode;
static int eno;
static float **rho,***rhox,***rhov;
static real dx[ND];
#ifdef PERIODIC
static int perbnd[ND];
#endif

/* gives the particle data (starting at offset 0) */


// Create trees for all particle types 0-5, use typelist[6] to control
// which particles to be combined if needed for tree generation
int treebuild_all(void)
{
    int i, tr;
    int nfree;
    int numnodes,startnode;

    startnode=1;
//    Box= All.BoxSize;
//    BoxHalf= All.BoxSize/2;
    nodes= nodes_base-startnode;
    snodes=snodes_base-startnode;
    if(All.TypeOfSmoothing==1)
	pnodes=pnodes_base;

    xnodes= xnodes_base-startnode;
    numnodestotal=0;
    pi=0;

#ifdef PERIODIC
    for(i=0;i<ND;i++)
    {
	if(All.boxh[i]>0)
	    perbnd[i]=1;
	if(All.boxh[i]<0)
	    perbnd[i]=-1;
	if(All.boxh[i]==0)
	    perbnd[i]=0;
    }
#endif

     trees.clear();
     for(i=0;i<int(npart.size());i++)
 	trees.push_back(0);


    // 'nfree' gives the next free internal tree node.
    for(tr=0, nfree=startnode, numnodestotal=0; tr<int(npart.size()); tr++)
    {
	if(npart[tr]>0)
	{
	    trees[tr]= nfree;
	    cout<<"\nParticle Type = "<<tr<<" First node number = "<<nfree<<endl;
	    nfree = treebuild_single(trees[tr], npart,tr, &numnodes);
	    numnodestree[tr]= numnodes;
	    cout<<"                  Last node number  = "<<nfree-1<<endl;
	}
	else
	{
	    trees[tr]= -1;
	    numnodestree[tr]= 0;
	}
    }
    for(i=1;i<nfree;++i)
	nodes[i].bnd=xnodes+i;


    cout<<"Total number of nodes = "<<numnodestotal<<endl;
    if(All.NodeSplittingCriterion==1)
    {
	if(All.CubicCells!=1)
	{
	    for(int i=0;i<ND;i++)
		delete [] rho[i];
	    delete [] rho;
	}

	if(All.CubicCells==1)
	{
	    //    int rows=int(pow(All.MaxPart,0.3333));
	    for(int i=0;i<2;i++)
	    {
		for(int j=0;j<2;j++)
		{
		    delete [] rhox[i][j];
		    delete [] rhov[i][j];
		}
		delete [] rhox[i];
		delete [] rhov[i];
	    }
	    delete [] rhox;
	    delete [] rhov;
	}
    }

    return numnodestotal;
}

// Create tree starting at index startnode with particles specified in typelist[6]
// creates a list of particles with given typelist
// and generates a root node for it
// output next available free node index, upadate total cumulative number of nodes created
// numnodestotal and provide nodes created for this typelist in creatednodes
int treebuild_single(int startnode, vector<int> &npart, int type,int *creatednodes)
{
    int  i, j,fp,fp1;
    int  nfree, numnodes;
    real x[ND][2],temp;
    struct NODE *nfreep;   struct SNODE *snfreep;  struct XNODE *xnfreep;
    nnodes= startnode;  nfreep= &nodes[nnodes];  snfreep= &snodes[nnodes];  xnfreep= &xnodes[nnodes];
    nfree=nnodes;
    /* find first and last particle of this type */
    for(i=0,fp=0; i<type; i++)
	fp+=npart[i];
    fp1=fp+npart[type];
    NumPart1=npart[type];

    rootnode=nnodes-1;
    list_kd[fp]=fp;
    nfreep->count=1;
    snfreep->lid=fp;
    nfreep->lid=fp;

#ifndef MEDIAN
    nfreep->parent=0;
#endif

    /* find enclosing rectangle */
    for(j=0; j<ND; j++)
	x[j][0]= x[j][1]= Part[fp].Pos[j];


    for(i=fp+1; i<fp1; i++)
	{
	    list_kd[i]=i;
	    nfreep->count++;
	    for(j=0; j<ND; j++)
	    {
		if(Part[i].Pos[j]>x[j][1])
		    x[j][1]= Part[i].Pos[j];
		if(Part[i].Pos[j]<x[j][0])
		    x[j][0]= Part[i].Pos[j];
	    }
	}


    /* determine maxmimum extension */
/*  if (All.PartBoundary==-1) temp=1e-6;
    else */
    //    temp=0.5e-6;
    //  temp=0.5/(pow(NumPart1*1.0,1.0/ND)+1);
    temp=0.5/(NumPart1*1.0-1.0);
    //  temp=2.0;
    if (All.PartBoundary==0) temp=1.0/(NumPart1*1.0-1.0);

    for(j=0; j<ND; j++)
    {
	kbuf[j]=(x[j][1]-x[j][0])*temp;
	xbmax[j]=x[j][1]+kbuf[j];
	xbmin[j]=x[j][0]-kbuf[j];
#ifdef PERIODIC
	if(perbnd[j]==0)
	{
	    if(All.boxh[j]<((xbmax[j]-xbmin[j])*0.5))
		All.boxh[j]=(xbmax[j]-xbmin[j])*0.5;
	}
	if(perbnd[j]==1)
	{
	    if((2*All.boxh[j])<(xbmax[j]-xbmin[j]-2*kbuf[j]))
	    {
		cout<<"Specified periodic length in "<<j<<" th dimension is too small"<<endl;
		cout<<"Try setting it to at least "<<(xbmax[j]-xbmin[j]-2*kbuf[j])*All.hs[j]<<endl;
		cout<<"Or let it be scaled automatically by setting it to 0.0"<<endl;
		endrun(10);
	    }
	}
	if(perbnd[j]==(-1))
	{
	    if(All.boxh[j]<((xbmax[j]-xbmin[j])*100))
		All.boxh[j]=(xbmax[j]-xbmin[j])*100;

	}
#endif

    }

    /* create a root node and insert first particle of correct type as its leaf */
    for(j=0; j<ND; j++)
    {
	xnfreep->x[j][1]=xbmax[j];
	xnfreep->x[j][0]=xbmin[j];
	snfreep->bleft[j]=1;
	snfreep->bright[j]=1;
    }


    snfreep->k0=ND;

    numnodes=numnodestotal;
    numnodestotal++;
    nfree1=nfree;
    treebuild(nfree);



    if(numnodestotal>=MaxNodes)
    {
	printf("maximum number %d of tree-nodes reached. \n", numnodestotal);
	printf("for particle %d\n", fp);
	endrun(1);
    }
    *creatednodes= numnodestotal-numnodes;

    i=1;j=1;
    while(i<nfreep->count)
    {
	i=i<<1;
	j=j+i;
    }

//  return nfree+numnodestotal-numnodes;
//  return nfree+j;
    return nfree1+1;
}


// Recursive routine to create a tree from an existing node
// needs strctures nodes,snodes,xnodees,pnodes and arrays list,list_kdt and father
void treebuild(int startnode)
{
    int i,j,k,l,ii,jl,jr,jmin,jmax,nfree;
    int pleft=0,pright=0,ks=0,k00,left,right;
    real kmean[ND],ksig[ND],kmax[ND],kmin[ND],m_tot;//,kmax1[ND],kmin1[ND];
    real kleft,kright,kcut,smin,smax;
    real ex=0,ev=0,e[ND];
    real temp=0.0,temp1=0.0,temp2=0.0,temp3=0.0;
    struct NODE *nfreep,*nfreepl,*nfreepr;   struct SNODE *snfreep,*snfreepl,*snfreepr;
    struct XNODE *xnfreep,*xnfreepl,*xnfreepr;
//  struct PNODE *pnfreep;
    nfree= startnode;
    nfreep= &nodes[nfree]; snfreep= &snodes[nfree]; xnfreep= &xnodes[nfree];
//  if(All.TypeOfSmoothing==1)
//  pnfreep= &pnodes[nfree];

//  nfree= startnode-rootnode;
//  cout<<startnode<<" "<<nfree<<endl;

    /* Calculate kmin.kmax,kmean */
    /*Check if leaf or not */
    if(nfreep->count > 1)
    {

#ifndef MEDIAN
	numnodestotal+=2;
#endif



	for(j=0; j<ND; j++)
	{
	    ksig[j] =0;
	    kmean[j]=0;
	}

	if(All.PartBoundary>1)
	{
	    i=list_kd[snfreep->lid];
	    for(j=0; j<ND; j++)
	    {
		kmin[j]=Part[i].Pos[j];
		kmax[j]=Part[i].Pos[j];
	    }
	}

	for(ii=0,m_tot=0;ii<nfreep->count;ii++)
	{
	    i=list_kd[snfreep->lid+ii];
	    for(j=0; j<ND; j++)
	    {
		kmean[j]+=Part[i].Pos[j];
	    }
	    m_tot+=Part[i].Mass;
	    //      if(All.PartBoundary>1)
	    for(j=0; j<ND; j++)
	    {
		if(kmin[j]>=Part[i].Pos[j])
		    kmin[j]=Part[i].Pos[j];
		if(kmax[j]<=Part[i].Pos[j])
		    kmax[j]=Part[i].Pos[j];
	    }
	}

	for(j=0; j<ND; j++)
	    kmean[j]=kmean[j]/nfreep->count;


// new boundary here
	/* Apply Boundary condition */


	if((All.PartBoundary>1)&&(nfreep->count>=All.PartBoundary))
	{
// decrease temp  to correct underestimation in  outer parts  
// decrease temp1 to correct underestimation in  inner parts   
	    if(All.CubicCells)
	    {
		temp=2.5*pow(NumPart1,1.0/ND);temp2=1.0*temp/5; temp3=2.0;
	    }
	    else
	    {
		temp=10*pow(NumPart1,1.0/ND); temp2=1.0*temp/5; temp3=2.0;

	    }

	    for(j=0;j<ND;j++)
	    {

#ifdef PERIODIC
	    if(perbnd[j]<0)
	    {
#endif

		if((kmax[j]-kmin[j])/(xnfreep->x[j][1]-xnfreep->x[j][0])> 1.0/NumPart1)
		{
		  
		    temp1=(kmax[j]-kmin[j])/(nfreep->count-1);
		    if(((xnfreep->x[j][1]-kmax[j])>(temp2*temp1))&&((kmin[j]-xnfreep->x[j][0])>(temp2*temp1)))
		    {
			xnfreep->x[j][1]=kmax[j]+temp3*temp1;
			xnfreep->x[j][0]=kmin[j]-temp3*temp1;
		    }
		    else
		    {
			if((xnfreep->x[j][1]-kmax[j])>(temp*temp1))
			    xnfreep->x[j][1]=kmax[j]+temp3*temp1;
			if((kmin[j]-xnfreep->x[j][0])>(temp*temp1))
			    xnfreep->x[j][0]=kmin[j]-temp3*temp1;
		    }

		}
#ifdef PERIODIC
	    }
#endif
	      

	    }
	}




/** Calculate the entropy */
	if((All.NodeSplittingCriterion==1)&&(All.CubicCells==1))
	    if(ND > 3)
	    {

		eno=int(pow(nfreep->count,0.333));
//        eno=int(pow(eno*1.0,0.333));
		if(nfreep->count<=8)
		    eno=2;

		for(i=0;i<eno;i++)
		    for(j=0;j<eno;j++)
			for(k=0;k<eno;k++)
			{
			    rhox[i][j][k]=0;
			    rhov[i][j][k]=0;
			}
		for(j=0;j<ND;j++)
		    dx[j]=(xnfreep->x[j][1]-xnfreep->x[j][0])/eno;
		//      cout<<nfree<<" "<<nfreep->count<<" "<<snfreep->ex<<" "<<snfreep->ev<<endl;
		for(ii=0;ii<nfreep->count;ii++)
		{
		    i=list_kd[snfreep->lid+ii];
		    j=int((Part[i].Pos[0]-xnfreep->x[0][0])/dx[0]);
		    k=int((Part[i].Pos[1]-xnfreep->x[1][0])/dx[1]);
		    l=int((Part[i].Pos[2]-xnfreep->x[2][0])/dx[2]);
		    rhox[j][k][l]+=Part[i].Mass;
		    j=int((Part[i].Pos[3]-xnfreep->x[3][0])/dx[3]);
		    k=int((Part[i].Pos[4]-xnfreep->x[4][0])/dx[4]);
		    l=int((Part[i].Pos[5]-xnfreep->x[5][0])/dx[5]);
		    rhov[j][k][l]+=Part[i].Mass;
		    //           cout<<j<<" "<<k<<" "<<l<<" rhov="<<rhov[j][k][l]<<" ii="<<ii<<" en="<<eno<<endl;
		}

		ex=0;ev=0;
		for(i=0;i<eno;i++)
		    for(j=0;j<eno;j++)
			for(k=0;k<eno;k++)
			{
			    if(rhox[i][j][k]>0)
			    {
				rhox[i][j][k]/=m_tot;
				ex-=rhox[i][j][k]*log10(rhox[i][j][k]);
			    }
			    if(rhov[i][j][k]>0)
			    {
				rhov[i][j][k]/=m_tot;
				ev-=rhov[i][j][k]*log10(rhov[i][j][k]);
			    }
			}
	    }


	/* choose the splitting dimension */
	if((All.NodeSplittingCriterion==1)&&(All.CubicCells!=1))
	{
	    eno=nfreep->count;
	    for(j=0;j<ND;j++)
		for(k=0;k<eno;k++)
		    rho[j][k]=0;

	    for(j=0;j<ND;j++)
		dx[j]=(xnfreep->x[j][1]-xnfreep->x[j][0])/eno;


	    for(ii=0;ii<nfreep->count;ii++)
	    {
		i=list_kd[snfreep->lid+ii];
		for(j=0;j<ND;j++)
		{
		    k=int((Part[i].Pos[j]-xnfreep->x[j][0])/dx[j]);
		    rho[j][k]+=Part[i].Mass;
		}
	    }
	    for(j=0;j<ND;j++)
		e[j]=0;
	    for(k=0;k<eno;k++)
	    {
		for(j=0;j<ND;j++)
		{
		    if(rho[j][k]>0)
		    {
			rho[j][k]/=m_tot;
			e[j]-=rho[j][k]*log10(rho[j][k]);
		    }
		}
	    }

	}




	for(ii=0;ii<nfreep->count;ii++)
	{
	    i=list_kd[snfreep->lid+ii];
	    for(j=0; j<ND; j++)
	    {
		ksig[j]+=(Part[i].Pos[j]-kmean[j])*(Part[i].Pos[j]-kmean[j]);
	    }

	}

	for(j=0; j<ND; j++)
	    ksig[j]=ksig[j]/nfreep->count;


//! check for spurious entropy of lattice like systems
	if((All.NodeSplittingCriterion==1)&&(All.CubicCells==1))
	{
	    jmin=0;      jmax=3;

	    if(ND >3)
	    {

		for(j=0,temp1=1.0;j<3;j++)
		    temp1*=((kmax[j]-kmin[j])/dx[j]);
		for(j=3,temp2=1.0;j<6;j++)
		    temp2*=((kmax[j]-kmin[j])/dx[j]);

		if(nfreep->count<=2)
		{temp1=2;temp2=2;}
		if((temp1<1)&&(temp2<1))
		    if(snfreep->k0 < 3){jmin=3;jmax=ND;}

		if((temp1>=1)&&(temp2>=1))
		{
		    if(ev==ex)
		    {
			if(snfreep->k0 < 3)
			{
			    jmin=3;                 jmax=ND;
			}
		    }
		    else
		    {
			if(ev<ex)
			{
			    jmin=3;
			    jmax=ND;
			}
		    }
		}
		if((temp2>=1)&&(temp1<1))
		{jmin=3;jmax=ND;}
	    }


	    for(j=jmin,smin=0.0,smax=ksig[jmin]; j<jmax; j++)
	    {
		if(ksig[j]>=smax)
		{
		    smax=ksig[j];                  ks=j;
		}
	    }
	}
	if((All.NodeSplittingCriterion==1)&&(All.CubicCells!=1))
	{
	    if(nfreep->count>=2)
	    {
		jmax=0;smin=-log10(1.0/nfreep->count);
		for(j=0; j<ND; j++)
		{
		    if((kmax[j]-kmin[j])/(xnfreep->x[j][1]-xnfreep->x[j][0])>=eno*0.01/NumPart1)
		    {
			if((e[j]==smin)&&(snfreep->k0==ks))  {smin=e[j]; ks=j;temp1=temp1+1;}
			if(e[j]<smin)  {smin=e[j]; ks=j;}
			jmax=1;
		    }
// 		    else
// 			cout<<" Lattice like structure or sparsely populated node, Check parameter PartBoundary "<<nfree<<" "<<j<<" "<<nfreep->count<<" "<<(kmax[j]-kmin[j])/(xnfreep->x[j][1]-xnfreep->x[j][0])<<endl;

		}

		if(jmax==0)
		{
		    cout<<" Warning:Lattice like structure or sparsely  populated node, Check parameter PartBoundary "<<nfree<<" "<<j<<" "<<nfreep->count<<" "<<smin<<endl;
		    k00=1+int(drand48()*ND);
		    while(k00==ND) k00=1+int(drand48()*ND);
		    ks=(snfreep->k0+k00)%ND;
		}
	    }
	    else
	    {
		ks=(snfreep->k0+1)%ND;
	    }
	}

	if(All.NodeSplittingCriterion==0)
	{
	    if(All.CubicCells==1)
	    {
		jmin=0;jmax=3;
		if((ND>3)&&(snfreep->k0 < 3))
		{ jmin=3;jmax=ND;}
		for(j=jmin,smin=0.0,smax=ksig[jmin]; j<jmax; j++)
		{
		    if(ksig[j]>=smax)
		    { smax=ksig[j];  ks=j; }
		}
	    }
	    if(All.CubicCells!=1)
	    {
		ks=(snfreep->k0+1)%ND;
		//    cout<<snfreep->k0<<" "<<ks<<" "<<(snfreep->k0+1)%ND-1<<endl;
	    }


	}

	// old boundary here


	int jj=0;
	jl=-1;
	jr=nfreep->count;
	while((jl==-1)||(jr==nfreep->count))
	{

	    k=ks;


	    /* Calculate kcut */
	    if(All.MedianSplittingOn==1)
	    {
		for(ii=0;ii<nfreep->count;ii++)
		{
		    i=list_kd[snfreep->lid+ii];
		    listr[snfreep->lid+ii]=Part[i].Pos[k];
		}
		float *lr;
		int * li;
		lr=listr+snfreep->lid-1;
		li=list_kd+snfreep->lid-1;
		sort2(nfreep->count,lr,li);
		jl=nfreep->count/2; jr=jl+1;
		kleft=*(lr+jl);   kright=*(lr+jr);
		jl--;jr--;
	    }
	    else
	    {
		kleft =xnfreep->x[k][0];
		kright=xnfreep->x[k][1];

		for(ii=0,jl=-1,jr=nfreep->count;ii<nfreep->count;ii++)
		{
		    i=list_kd[snfreep->lid+ii];
		    if(Part[i].Pos[k]<= kmean[k])
		    {
			if(Part[i].Pos[k] >= kleft)
			{
			    kleft=Part[i].Pos[k];
			    pleft=i;
			}
			list_kdt[++jl]=i;
		    }
		    else
		    {
			if(Part[i].Pos[k]<= kright)
			{
			    kright=Part[i].Pos[k];
			    pright=i;
			}
			list_kdt[--jr]=i;
		    }
		}


	    }
	    //------------------------------------------------------


	    if(jl==-1)
	    {
//corrected  2007 Jan 11
//		i=pright;
		i=list_kdt[jr];
		cout<<"particles too close left par no = 0  "<<k<<" "<<numnodestotal<<" "<<kmean[k]<<" "<<kleft<<" "<<nfreep->count<<endl;
		if(jj==ND)
		{
#ifdef WARN
		    checkrun("Particles have identical co-ordinates","exit","correct it");
#else
		    cout<<"Particles have identical co-ordinates correcting it"<<endl;
#endif
// bug fixed
		    Part[i].Pos[k]=xnfreep->x[k][0]+(kmean[k]-xnfreep->x[k][0])/2.0;
		    kleft=Part[i].Pos[k];
		    jl++;
		    jr++;
		}
		else
		    ks=(ks+1)%ND;

	    }

	    if(jr==nfreep->count)
	    {
//corrected 2007 Jan 11
//		i=pleft;
		i=list_kdt[jl];

		cout<<"particles too close right par no = 0  "<<k<<" "<<numnodestotal<<" "<<kmean[k]<<" "<<kleft<<" "<<nfreep->count<<endl;
		if(jj==ND)
		{
#ifdef WARN
		    checkrun("WARNING! Particles have identical co-ordinates","exit","correct it");
#else
		    cout<<"Particles have identical co-ordinates correcting it"<<endl;
#endif
// bug fixed
		    Part[i].Pos[k]=xnfreep->x[k][1]-(xnfreep->x[k][1]-kmean[k])/2.0;
		    kright=Part[i].Pos[k];
		    jr--;
		    jl--;
		}
		else
		    ks=(ks+1)%ND;

	    }

	    jj++;
	}


	if(All.MedianSplittingOn!=1)
	    for(ii=0;ii<nfreep->count;ii++)
		list_kd[snfreep->lid+ii]=list_kdt[ii];



	kcut=(kleft+kright)/2.0;

	//      cout<<nfree<<" "<<ks<<" "<<kcut<<nfreep->count<<" "<<ev<<" "<<ex<<endl;

	/* create and update the left right nodes */
	nfreep->kcut=kcut;
	nfreep->k1=k;





#ifdef MEDIAN
 	left= LOWER(nfree,nodes1)+rootnode;
 	right = UPPER(nfree,nodes1)+rootnode;
#else
	++nnodes;left= nnodes;
	nfreep->left= left-rootnode;
	++nnodes;right = nnodes;
#endif
	nfreepl=&nodes[left];    snfreepl=&snodes[left];  xnfreepl=&xnodes[left];
	nfreepr=&nodes[right];   snfreepr=&snodes[right]; xnfreepr=&xnodes[right];


	if(nfree1<right) nfree1=right;
#ifdef MEDIAN
	numnodestotal=nfree1;
#endif

	for(j=0;j<ND;j++)
	{
	    if(j==k)
	    {
		xnfreepl->x[j][1]=kcut;
		xnfreepl->x[j][0]=xnfreep->x[j][0];
		xnfreepr->x[j][1]=xnfreep->x[j][1];
		xnfreepr->x[j][0]=kcut;
		snfreepl->bright[j]=0;
		snfreepl->bleft[j]=snfreep->bleft[j];
		snfreepr->bleft[j]=0;
		snfreepr->bright[j]=snfreep->bright[j];
	    }
	    else
	    {
		xnfreepl->x[j][1]=xnfreep->x[j][1];
		xnfreepl->x[j][0]=xnfreep->x[j][0];
		xnfreepr->x[j][1]=xnfreep->x[j][1];
		xnfreepr->x[j][0]=xnfreep->x[j][0];
		snfreepl->bleft[j] =snfreep->bleft[j];
		snfreepl->bright[j]=snfreep->bright[j];
		snfreepr->bright[j]=snfreep->bright[j];
		snfreepr->bleft[j] =snfreep->bleft[j];
	    }
	}
	snfreepl->k0=k;
	snfreepr->k0=k;


#ifndef MEDIAN
	{
	    nfreepl->left=NULL_N;
	    nfreepr->left=NULL_N;
	    nfreepl->parent=nfree-rootnode;
	    nfreepr->parent=nfree-rootnode;
	}
#endif
	snfreepl->lid=snfreep->lid;
	snfreepr->lid=snfreep->lid+jl+1;

	nfreepl->lid=snfreep->lid;
	nfreepr->lid=snfreep->lid+jl+1;

	nfreepl->count=jl+1;
	nfreepr->count=nfreep->count-jl-1;
	treebuild(left);
	treebuild(right);
    }
    else
    {
	/* if leaf wind up */
	nfreep->kcut=0.0;
	i=list_kd[snfreep->lid];

	if(All.PartBoundary==1)
	    for(j=0; j<ND; j++)
	    {
#ifdef PERIODIC
	    if(perbnd[j]<0)
	    {
#endif
		if(snfreep->bright[j]==1)
		    xnfreep->x[j][1]=Part[i].Pos[j]+(Part[i].Pos[j]-xnfreep->x[j][0]);
		if(snfreep->bleft[j]==1)
		    xnfreep->x[j][0]=Part[i].Pos[j]-(xnfreep->x[j][1]-Part[i].Pos[j]);
#ifdef PERIODIC
	    }
#endif
	    }

	if(All.TypeOfSmoothing==1)
	    for(j=0;j<ND;j++)
	    {
		pnodes[pi].x[j][0]=xnfreep->x[j][0]; pnodes[pi].x[j][1]=xnfreep->x[j][1];
	    }
	nfreep->k1=pi;        pi++;

    }
}


// Allocate memory
void treeallocate(int maxpart)  /* for median maxnodes=4*maxpart is sufficient */
{
    int bytes=0,i,j,rows;
    float allbytes=0;
#ifdef MEDIAN
    MaxNodes=4*maxpart;
#else
    MaxNodes=2*maxpart;
#endif

    if(!(list_kd=new int4byte[maxpart]))
    {
	fprintf(stdout,"Failed to allocate %d spaces for 'list' array (%d bytes)\n", maxpart, bytes);
	exit(0);
    }
    allbytes+=maxpart*sizeof(int4byte);

    if(!(list_kdt=new int4byte[maxpart]))
    {
	fprintf(stdout,"Failed to allocate %d spaces for 'list_kdt' array (%d bytes)\n", maxpart, bytes);
	exit(0);
    }
    allbytes+=maxpart*sizeof(int4byte);


    if(!(nodes_base=new NODE[MaxNodes+1]))
    {
	printf("failed to allocate memory for %d tree-nodes (%d bytes).\n", MaxNodes, bytes);
	endrun(3);
    }
    allbytes+=MaxNodes*sizeof(struct NODE);


    if(!(xnodes_base=new XNODE[MaxNodes+1]))
    {
	printf("failed to allocate memory for %d tree-nodes (%d bytes).\n", MaxNodes, bytes);
	endrun(3);
    }
    allbytes+=MaxNodes*sizeof(struct XNODE);


    if(!(snodes_base=new SNODE[MaxNodes+1]))
    {
	printf("failed to allocate memory for %d tree-nodes (%d bytes).\n", MaxNodes, bytes);
	endrun(3);
    }
    allbytes+=MaxNodes*sizeof(struct SNODE);


    if(All.TypeOfSmoothing==1)
    {
	if(!(pnodes_base=new PNODE[maxpart+1]))
	{
	    printf("failed to allocate memory for %d tree-nodes (%d bytes).\n", MaxNodes, bytes);
	    endrun(3);
	}
	allbytes+=maxpart*sizeof(struct PNODE);
    }

    if(All.NodeSplittingCriterion==1)
    {
	if(All.CubicCells!=1)
	{
	    rho=new float* [ND];
	    for(i=0;i<ND;i++)
		if(!(rho[i]=new float [maxpart]))
		{
		    fprintf(stdout,"Failed to allocate %d spaces for 'entropy' bins (%d bytes)\n", maxpart, bytes);
		    exit(0);
		}
	    allbytes+=maxpart*ND*sizeof(float);
	}

	if(All.CubicCells==1)
	{
	    rows=int(pow(maxpart,0.3333));
	    rhox=new float** [rows];
	    for(i=0;i<rows;i++)
		if(!(rhox[i]=new float* [rows]))
		{
		    fprintf(stdout,"Failed to allocate %d spaces for 'entropy' bins (%d bytes)\n", maxpart, bytes);
		    exit(0);
		}
		else
		{
		    for(j=0;j<rows;j++)
			if(!(rhox[i][j]=new float [rows]))
			{
			    fprintf(stdout,"Failed to allocate %d spaces for 'entropy' bins (%d bytes)\n", maxpart, bytes);
			    exit(0);
			}
		}
	    allbytes+=rows*rows*rows*sizeof(float);

	    rhov=new float** [rows];
	    for(i=0;i<rows;i++)
		if(!(rhov[i]=new float* [rows]))
		{
		    fprintf(stdout,"Failed to allocate %d spaces for 'entropy' bins (%d bytes)\n", maxpart, bytes);
		    exit(0);
		}
		else
		{
		    for(j=0;j<rows;j++)
			if(!(rhov[i][j]=new float [rows]))
			{
			    fprintf(stdout,"Failed to allocate %d spaces for 'entropy' bins (%d bytes)\n", maxpart, bytes);
			    exit(0);
			}
		}
	    allbytes+=rows*rows*rows*sizeof(float);
	}
    }

    if(!(listr=new float[maxpart]))
    {
	fprintf(stdout,"Failed to allocate %d spaces for 'list' array (%d bytes)\n", maxpart, bytes);
	exit(0);
    }
    allbytes+=maxpart*sizeof(float);

//    cout<<" a "<<maxpart<<MaxNodes<<endl;




    printf("Allocated %g MByte for Binary-Tree and Entropy bins.\n\n",allbytes/(1024.0*1024.0));
    //  cout<<MaxNodes<<" "<<maxpart<<endl;

}


// free the allocated memory
void treefree(void)
{
    delete [] nodes_base;
    delete [] snodes_base;
    delete [] xnodes_base;
    delete [] list_kd;
    delete [] list_kdt;
    delete [] listr;
    delete [] pnodes_base;
}



