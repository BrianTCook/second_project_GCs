/***************************************************************************
                          tree_search.cpp  -  description
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
//#include "proto.h"
#include "nr.h"
#include "tree.h"
#include "functions.h"

static real msmooth,vsmooth;


// **************************************
//  *    Below we have the routines       *
//  *    dealing with neighbour finding.  *
//  **************************************


//  This routine maps a coordinate difference
//  to the nearest periodic image
//
void ngb_calculate_metric(float xyz[ND], struct NODE* nodes1)
{
    real sr[ND],fac=1.2,mp1;  /* search radius */
    int   th=0,j,it=0;

    nodesC=nodes1;
    mp1=1;
    if(searchcenter[0]==0.0)
    {
	th=1;
	while(nodes1[th].count>1)
	{
	    j=nodes1[th].k1;
	    if(xyz[j] < nodes1[LOWER(th,nodes1)].bnd->x[j][1])
		th=LOWER(th,nodes1);
	    else
		th=UPPER(th,nodes1);
	}

	for(j=0; j<ND; j++)
	{
	    sr[j]=real(nodes1[th].bnd->x[j][1]-nodes1[th].bnd->x[j][0])*0.5;
	    searchcenter[j]=sr[j];
	}

    }
    else
	for(j=0; j<ND; j++)
	{
	    sr[j]=searchcenter[j];
	    searchcenter[j]=0.0;
	}



    //  temp=(mp1+mp2)*0.5;
    for(j=0; j<ND; j++)
    {
	sr[j]=fac*sr[j]*0.5;
	//            sr[j]*=All.SearchFactor*0.5*pow(temp,1.0/ND);
	searchx[j][0]=xyz[j];
	searchx[j][1]=xyz[j];
    }
//    cout<<th<<" a "<<msmooth<<" "<<mp1<<" "<<" "<<sr[0]<<endl;
    msmooth=0.0; vsmooth=0.0;
//    msmooth=2;
    while(msmooth<mp1)
    {
	msmooth=0.0; vsmooth=0.0;
	for(j=0; j<ND; j++)
	{
	    sr[j]*=2;
	    searchx[j][0]=xyz[j]-sr[j];
	    searchx[j][1]=xyz[j]+sr[j];
	    searchcenter[j]=0.0;
	}

	ngb_treesearch_kdts_h(1);

// 	if(msmooth<0.1)
// 	{
// 	    msmooth=mp1;
// 		for(j=0; j<ND; j++)
// 		    searchcenter[j]=1.0;
// 		endrun(10);
// 	}
// 	else
	if(msmooth>mp1)
	{
		for(j=0; j<ND; j++)
		{
		    searchcenter[j]/=vsmooth;
		}


	}
	it++;
	if(it>100)
	{
	    cout<<it<<" too many iterations ngb_calculate _metric"<<msmooth<<endl;
	    for(j=0;j<ND;j++)
		cout<<xyz[j]<<"    "<<sr[j]<<endl;
	    endrun(2);
	}
	if(it>20)
	    cout<<th<<" a "<<msmooth<<" "<<mp1<<" "<<xyz[0]<<" "<<sr[0]<<endl;

//	    cout<<th<<" a "<<msmooth<<" "<<mp1<<" "<<xyz[0]<<" "<<sr[0]<<endl;

    }

}

// Fiestas style density calculation. Calculate a search box
//enclosing a given mass and then calculate density

float ngb_treedensity_fiestas(float xyz[ND], struct NODE* nodes1)
{
    double temp,temp1;
    real sr[ND],fac=1.0;  /* search radius */
    int   th=0,it,j;
    double mp1,mp2,mass_unit;

//    cout<<xyz[0]<<endl;
    mp1=All.DesNumNgb*0.99;
    mp2=All.DesNumNgb*1.01;


//    cout<<list_kd[nodes1[1].lid]<<endl;
    mass_unit=Part[list_kd[nodes1[1].lid]].Mass;



//    nodes1=nodes+(parttype-1);
    nodesC=nodes1;

    if(All.CubicCells==1)
    {
	th=1;
	while(nodes1[th].count>1)
	{
	    j=nodes1[th].k1;
	    if(xyz[j] < nodes1[LOWER(th,nodes1)].bnd->x[j][1])
		th=LOWER(th,nodes1);
	    else
		th=UPPER(th,nodes1);

	}

	for(j=0; j<ND; j++)
	{
	    sr[j]=(nodes1[th].bnd->x[j][1]-nodes1[th].bnd->x[j][0])*0.1*pow(double(All.DesNumNgb*1.0/nodes1[th].count),1.0/ ND);
	}


	for(j=0,temp=1.0; j<3; j++)
	    temp*=(sr[j]);
	temp=pow(temp,0.3333);
	for(j=0; j<3; j++)
	    sr[j]=temp;

	if(ND>3)
	{
	    for(j=3,temp=1.0; j<6; j++)
		temp*=(sr[j]);
	    temp=pow(temp,0.3333);
	    for(j=3; j<6; j++)
		sr[j]=temp;
	}
    }
    else
    {
	for(j=0; j<ND; j++)
	    searchcenter[j]=0.0;

	ngb_calculate_metric(xyz,nodes1);
	for(j=0; j<ND; j++)
	{
	    sr[j]=searchcenter[j]*0.1*pow(double(All.DesNumNgb*1.0/1.0),1.0/ ND);
	}

    }


    //  temp=(mp1+mp2)*0.5;
    for(j=0; j<ND; j++)
    {
	sr[j]=fac*sr[j]/2.0;
	//            sr[j]*=All.SearchFactor*0.5*pow(temp,1.0/ND);
	searchx[j][0]=xyz[j];
	searchx[j][1]=xyz[j];
    }
    //   cout<<th<<" a "<<msmooth<<" "<<mp1<<" "<<mp2<<" "<<msmooth<<endl;
    it=0;
    msmooth=0.0;
    vsmooth=0.0;
    while(msmooth<mp1)
    {
	//   cout<<th<<" b "<<msmooth<<" "<<mp1<<" "<<mp2<<" "<<sr[0]<<" "<<par<<endl;
	it++;
	if(it > 100)
	{
	    cout<<"Too many iterations check smoothing"<<endl;
	    cout<<th<<" "<<msmooth<<" "<<mp1<<" "<<mp2<<" "<<endl;
	    endrun(10);
	}
	msmooth=0.0;
	vsmooth=0.0;
	for(j=0; j<ND; j++)
	{
	    sr[j]*=2;
	    searchx[j][0]=xyz[j]-sr[j];
	    searchx[j][1]=xyz[j]+sr[j];
	}

	if(All.VolCorr==1)
	    ngb_treesearch_kdts_vc(1);
	else
	    ngb_treesearch_kdts(1);

	if((it==1)&&(msmooth>mp1))
	{
	    for(j=0; j<ND; j++)
		sr[j]/=4; // actually this halfs the box length since it was already doubled earlier
	    it=0; msmooth=0.0;
//	          cout<<th<<endl;
	}

    }

    for(j=0; j<ND; j++)
	sr[j]=sr[j]*0.5;
    //    cout<<th<<" c "<<msmooth<<" "<<mp1<<" "<<mp2<<" "<<msmooth<<endl;

    it=0;
    while((msmooth<mp1)||(msmooth>mp2))
    {
	//         cout<<th<<"  d "<<trees[parttype]<<" "<<sr[0]<<" "<<mp1<<" "<<mp2<<" "<<msmooth<<endl;

	it++;
	if(it > 100)
	{
	    cout<<"Too many iterations check smoothing"<<endl;
	    cout<<th<<" "<<msmooth<<" "<<mp1<<" "<<mp2<<" "<<endl;
	    endrun(10);
	}

	if(msmooth>mp2)
	    for(j=0; j<ND; j++)
	    {
		sr[j]=sr[j]*0.5;
		searchx[j][0]+=sr[j];
		searchx[j][1]-=sr[j];
	    }
	if(msmooth<mp1)
	    for(j=0; j<ND; j++)
	    {
		sr[j]=sr[j]*0.5;
		searchx[j][0]-=sr[j];
		searchx[j][1]+=sr[j];
	    }
	msmooth=0.0;
	vsmooth=0.0;

	if(All.VolCorr==1)
	    ngb_treesearch_kdts_vc(1);
	else
	    ngb_treesearch_kdts(1);

    }



    for(j=0,temp=1.0; j<ND; j++)
	temp*=(searchx[j][1]-searchx[j][0]);

    if(All.CubicCells==1) temp1=0.5; else temp1=0.9; //0.2 0.5

    if(All.VolCorr==1)
    {
	if (vsmooth/temp <temp1)
	    return (msmooth*mass_unit/vsmooth);
	else
	    return (msmooth*mass_unit/temp);
    }
    else
	return (msmooth*mass_unit/temp);

}


float ngb_treedensity_raw(float xyz[ND], struct NODE* nodes1)
{
    double temp;
    int   th,j;
//    nodes1=nodes+parttype-1;
    th=1;
    while(nodes1[th].count>1)
    {
	j=nodes1[th].k1;
	if(xyz[j] < nodes1[LOWER(th,nodes1)].bnd->x[j][1])
	    th=LOWER(th,nodes1);
	else
	    th=UPPER(th,nodes1);
    }

    for(j=0,temp=1.0; j<ND; j++)
    {
	temp*=real(nodes1[th].bnd->x[j][1]-nodes1[th].bnd->x[j][0]);
    }

    return Part[list_kd[nodes1[th].lid]].Mass/temp;

}


// Search kd tree to find particles in a given box
void ngb_treesearch_kdts(int no)
{
    int k,l;
    real kvol=1.0,temp1,temp2,temp3,temp4;
    struct NODE *nfreep;
    struct PNODE *pnfreep;
//    struct XNODE *pnfreep;

    nfreep= &nodesC[no];

    if(nfreep->count > 1)
    {
	k=nfreep->k1;
	temp1=nfreep->kcut;
	l=LOWER(no,nodesC);
	if(temp1>searchx[k][0])
	    ngb_treesearch_kdts(l);
	l=UPPER(no,nodesC);
	if(temp1<searchx[k][1])
	    ngb_treesearch_kdts(l);
    }
    else
    {
	pnfreep= &pnodes[nfreep->k1];
//	pnfreep=nodesC[no].bnd;
//	pnfreep= &xnodes[no];
	for(k=0; k<ND; k++)
	{
	    temp1=searchx[k][0];    temp2=searchx[k][1];
	    temp3=pnfreep->x[k][0]; temp4=pnfreep->x[k][1];
	    if(temp3>temp1)    temp1=temp3;
	    if(temp4<temp2)    temp2=temp4;
	    if(temp2<=temp1)   return;
	    else
	    {
		kvol*=(temp2-temp1)/(temp4-temp3);
	    }
	}
	msmooth+=kvol;
    }

    return;

}


// Search kd tree to find particles in a given box along with volume calculation
void ngb_treesearch_kdts_vc(int no)
{
    int k,l;
    real kvol=1.0,temp1,temp2,temp3,temp4;
    real kvol1=1.0;
    struct NODE  *nfreep;
    struct PNODE *pnfreep;
//    struct XNODE *pnfreep;
    nfreep= &nodesC[no];

    if(nfreep->count > 1)
    {
	k=nfreep->k1;
	temp1=nfreep->kcut;
	l=LOWER(no,nodesC);
	if(temp1>searchx[k][0])
	    ngb_treesearch_kdts_vc(l);
	l=UPPER(no,nodesC);
	if(temp1<searchx[k][1])
	    ngb_treesearch_kdts_vc(l);
    }
    else
    {
	pnfreep= &pnodes[nfreep->k1];
//	pnfreep=nodesC[no].bnd;
//	pnfreep= &xnodes[no];
	for(k=0; k<ND; k++)
	{
	    temp1=searchx[k][0];    temp2=searchx[k][1];
	    temp3=pnfreep->x[k][0]; temp4=pnfreep->x[k][1];
	    if(temp3>temp1)    temp1=temp3;
	    if(temp4<temp2)    temp2=temp4;
	    if(temp2<=temp1)   return;
	    else
	    {
		kvol*=(temp2-temp1)/(temp4-temp3);
		kvol1*=(temp2-temp1);
	    }
	}

	msmooth+=kvol;
	vsmooth+=kvol1;
    }
    return;
}



void ngb_treesearch_kdts_h(int no)
{
    int k,l;
    real temp1,temp2,temp3,temp4,kvol=1.0;
    struct NODE  *nfreep;
//    struct PNODE *pnfreep;
    struct XNODE *pnfreep;


    nfreep= &nodesC[no];

    if(nfreep->count > 1)
    {
	k=nfreep->k1;
	temp1=nfreep->kcut;

	l=LOWER(no,nodesC);
       	if(temp1>searchx[k][0])
	    ngb_treesearch_kdts_h(l);
	l=UPPER(no,nodesC);	
	if(temp1<searchx[k][1])
	    ngb_treesearch_kdts_h(l);
    }
    else
    {
//	pnfreep= &pnodes[nfreep->k1];
	pnfreep=nodesC[no].bnd;
//	pnfreep= &xnodes[no];
	for(k=0; k<ND; k++)
	{
	    temp1=searchx[k][0];    temp2=searchx[k][1];
	    temp3=pnfreep->x[k][0]; temp4=pnfreep->x[k][1];
	    if(temp3>temp1)    temp1=temp3;
	    if(temp4<temp2)    temp2=temp4;
	    if(temp2<=temp1)   return;
	    else
	    {
		kvol*=(temp2-temp1)/(temp4-temp3);
//		kvol1*=(temp2-temp1);
	    }
	}

	for(k=0; k<ND; k++)
	{
	    temp3=pnfreep->x[k][0]; temp4=pnfreep->x[k][1];
	    searchcenter[k]+=((temp4-temp3));
//	    cout<<k<<" "<<searchcenter[k]<<endl;
	}
	msmooth+=kvol;
	vsmooth=vsmooth+1.0;
    }
    return;
}



void treeprint(int no)
{
    int k,l;
    //  real kvol=1.0,temp1,temp2,temp3,temp4;
    struct NODE  *nfreep;
    //  struct PNODE *pnfreep;
    struct SNODE *snfreep;
    struct XNODE *xnfreep;

    nfreep= &nodes[no];
    xnfreep= &xnodes[no];
    snfreep= &snodes[no];
    fprintf(treefile_fpt," %d %e %d %d\n",nfreep->k1,nfreep->kcut,nfreep->count,snfreep->lid);
    for(k=0;k<ND;k++)
	fprintf(treefile_fpt," %e %e\n",xnfreep->x[k][0],xnfreep->x[k][1]);

    if(nfreep->count > 1)
    {
	l=LOWER(no,nodesC);
	treeprint(l);
	l=UPPER(no,nodesC);
	treeprint(l);
    }
    return;
}




#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define SWAPI(a,b) tempi=(a);(a)=(b);(b)=tempi;
#define M 7
#define NSTACK 50
void sort2(unsigned long n, float arr[], int brr[])
{
    unsigned long i,ir=n,j,k,l=1;
    int *istack,jstack=0;
    float a,temp;
    int b,tempi;
    istack=ivector(1,NSTACK);
    for (;;) {
	if (ir-l < M) {
	    for (j=l+1;j<=ir;j++) {
		a=arr[j];
		b=brr[j];
		for (i=j-1;i>=l;i--) {
		    if (arr[i] <= a) break;
		    arr[i+1]=arr[i];
		    brr[i+1]=brr[i];
		}
		arr[i+1]=a;
		brr[i+1]=b;
	    }
	    if (!jstack) {
		free_ivector(istack,1,NSTACK);
		return;
	    }
	    ir=istack[jstack];
	    l=istack[jstack-1];
	    jstack -= 2;
	} else {
	    k=(l+ir) >> 1;
	    SWAP(arr[k],arr[l+1])
		SWAPI(brr[k],brr[l+1])
		if (arr[l] > arr[ir]) {
		    SWAP(arr[l],arr[ir])
			SWAPI(brr[l],brr[ir])
			}
	    if (arr[l+1] > arr[ir]) {
		SWAP(arr[l+1],arr[ir])
		    SWAPI(brr[l+1],brr[ir])
		    }
	    if (arr[l] > arr[l+1]) {
		SWAP(arr[l],arr[l+1])
		    SWAPI(brr[l],brr[l+1])
		    }
	    i=l+1;
	    j=ir;
	    a=arr[l+1];
	    b=brr[l+1];
	    for (;;) {
		do i++; while (arr[i] < a);
		do j--; while (arr[j] > a);
		if (j < i) break;
		SWAP(arr[i],arr[j])
		    SWAPI(brr[i],brr[j])
		    }
	    arr[l+1]=arr[j];
	    arr[j]=a;
	    brr[l+1]=brr[j];
	    brr[j]=b;
	    jstack += 2;
	    if (jstack > NSTACK) nrerror("NSTACK too small in sort2.");
	    if (ir-i+1 >= j-l) {
		istack[jstack]=ir;
		istack[jstack-1]=i;
		ir=j-1;
	    } else {
		istack[jstack]=j-1;
		istack[jstack-1]=l;
		l=i;
	    }
	}
    }
}
#undef M
#undef NSTACK
#undef SWAP
#undef SWAPI


